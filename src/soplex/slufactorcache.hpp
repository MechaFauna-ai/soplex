// SPDX-License-Identifier: Apache-2.0

#ifndef _SLUFACTORCACHE_HPP_
#define _SLUFACTORCACHE_HPP_

#include <algorithm>
#include <cmath>
#include <new>
#include <utility>

namespace soplex
{
template <class R>
bool SLUFactor<R>::factorCacheEnabled() const
{
#if defined(__SSE2__) && FLT_EVAL_METHOD == 0
   return std::is_same<R, double>::value && std::is_same<Real, double>::value
      && this->thedim >= 128 && this->thedim <= 65536
      && std::fegetround() == FE_TONEAREST && (_mm_getcsr() & 0xffc0) == 0x1f80;
#else
   return false;
#endif
}

template <class R>
void SLUFactor<R>::factorCacheClear()
{
   for(auto& slot : factorSlots) slot = FactorSlot();
   factorTick = 0;
}

template <class R>
int SLUFactor<R>::factorCacheBegin(const SVectorBase<R>* matrix[], int dm)
{
   if(!factorCacheEnabled()) return -1;
   // These addresses only screen likely repeats. No stale pointer is dereferenced.
   uint64_t hash = uint64_t(dm);
   for(int j = 0; j < dm; ++j)
      hash = (hash ^ uint64_t(reinterpret_cast<uintptr_t>(matrix[j]))) * UINT64_C(1099511628211);
   int oldest = 0;
   for(int i = 0; i < 4; ++i)
   {
      auto& slot = factorSlots[i];
      if(slot.hash == hash && slot.pointers.size() == size_t(dm)
         && std::equal(slot.pointers.begin(), slot.pointers.end(), matrix))
      {
         slot.tick = ++factorTick;
         if(slot.visits < 2) ++slot.visits;
         return i;
      }
      if(slot.tick < factorSlots[oldest].tick) oldest = i;
   }
   try
   {
      auto& slot = factorSlots[oldest];
      slot = FactorSlot();
      slot.pointers.assign(matrix, matrix + dm);
      slot.hash = hash;
      slot.tick = ++factorTick;
      slot.visits = 1;
   }
   catch(const std::bad_alloc&) { return -1; }
   return oldest;
}

template <class R>
bool SLUFactor<R>::factorCacheKey(const SVectorBase<R>* matrix[], std::vector<uint64_t>& key) const
{
   // Require the full ordered input and every numerical LU control, including
   // signed-zero coefficient bits. A hash or pointer match is never sufficient.
   const auto bits = [](double value) { uint64_t out; std::memcpy(&out, &value, sizeof(out)); return out; };
   key.clear();
   size_t words = 12 + size_t(this->thedim);
   for(int j = 0; j < this->thedim; ++j)
   {
      words += size_t(matrix[j]->size()) * 2;
      if(words > 2 * 1024 * 1024) return false;
   }
   key.reserve(words);
   key.push_back(this->thedim);
   key.push_back(bits(double(lastThreshold)));
   key.push_back(bits(double(this->tolerances()->epsilonFactorization())));
   key.push_back(bits(double(this->tolerances()->epsilonPivot())));
   key.push_back(bits(double(this->tolerances()->epsilonUpdate())));
   key.push_back(bits(double(this->tolerances()->epsilon())));
   key.push_back(this->l.updateType);
   key.push_back(bits(double(this->rowMemMult)));
   key.push_back(bits(double(this->colMemMult)));
   key.push_back(bits(double(this->lMemMult)));
   for(int j = 0; j < this->thedim; ++j)
   {
      const auto& col = *matrix[j];
      key.push_back(col.size());
      for(int k = 0; k < col.size(); ++k)
      {
         double value = double(col.value(k));
         if(!std::isfinite(value)) return false;
         key.push_back(col.index(k));
         key.push_back(bits(value));
      }
   }
   return true;
}

template <class R>
std::shared_ptr<typename SLUFactor<R>::FactorImage> SLUFactor<R>::factorImage() const
{
   // Count allocated array lengths as well as values, including update capacity.
   const size_t bytes = (this->diag.size() + this->u.row.val.size() + this->u.col.val.size()
      + this->l.val.size() + this->l.rval.size()) * sizeof(R)
      + (size_t(this->u.row.size) + this->u.col.size + this->l.size
         + size_t(this->thedim) * 16 + size_t(this->l.startSize) * 2
         + this->l.rval.size()) * sizeof(int);
   // Bound estimated storage per image; at most two images remain admitted.
   if(bytes > 64 * 1024 * 1024) return nullptr;
   if(!std::isfinite(double(this->initMaxabs)) || !std::isfinite(double(this->maxabs))) return nullptr;
   for(const auto& value : this->diag) if(!std::isfinite(double(value))) return nullptr;
   auto image = std::make_shared<FactorImage>();
   auto& x = *image;
   x.diag = this->diag; x.urval = this->u.row.val; x.ucval = this->u.col.val;
   x.lval = this->l.val; x.lrval = this->l.rval;
   const auto copy = [](std::vector<int>& out, const int* in, size_t n)
   {
      out.resize(n);
      if(n) std::memcpy(out.data(), in, n * sizeof(int));
   };
   const int n = this->thedim;
   copy(x.rp, this->row.perm, n); copy(x.ro, this->row.orig, n);
   copy(x.cp, this->col.perm, n); copy(x.co, this->col.orig, n);
   copy(x.uridx, this->u.row.idx, this->u.row.size);
   copy(x.urstart, this->u.row.start, n + 1); copy(x.urlen, this->u.row.len, n + 1);
   copy(x.urmax, this->u.row.max, n + 1);
   copy(x.ucidx, this->u.col.idx, this->u.col.size);
   copy(x.ucstart, this->u.col.start, n + 1); copy(x.uclen, this->u.col.len, n + 1);
   copy(x.ucmax, this->u.col.max, n + 1);
   for(auto* p = this->u.row.list.next; p != &this->u.row.list; p = p->next) x.urorder.push_back(p->idx);
   for(auto* p = this->u.col.list.next; p != &this->u.col.list; p = p->next) x.ucorder.push_back(p->idx);
   copy(x.lidx, this->l.idx, this->l.size);
   copy(x.lstart, this->l.start, this->l.startSize); copy(x.lrow, this->l.row, this->l.startSize);
   if(this->l.rbeg)
   {
      copy(x.lridx, this->l.ridx, this->l.start[this->l.firstUpdate]);
      copy(x.lrbeg, this->l.rbeg, n + 1);
      copy(x.lrorig, this->l.rorig, n); copy(x.lrperm, this->l.rperm, n);
   }
   x.urused = this->u.row.used; x.ucused = this->u.col.used;
   x.firstUpdate = this->l.firstUpdate; x.firstUnused = this->l.firstUnused;
   x.updateType = this->l.updateType; x.nz = this->nzCnt;
   x.initMaxabs = this->initMaxabs; x.maxabs = this->maxabs;
   return image;
}

template <class R>
void SLUFactor<R>::restoreFactorImage(const FactorImage& x)
{
   const auto alloc = [](int*& out, const std::vector<int>& in, size_t oldSize = 0)
   {
      if(in.size() > oldSize) spx_realloc(out, in.size());
   };
   // Complete every potentially throwing allocation before replacing contents.
   this->diag.reserve(x.diag.size()); this->u.row.val.reserve(x.urval.size());
   this->u.col.val.reserve(x.ucval.size()); this->l.val.reserve(x.lval.size());
   this->l.rval.reserve(x.lrval.size());
   alloc(this->u.row.idx, x.uridx, this->u.row.size); alloc(this->u.col.idx, x.ucidx, this->u.col.size);
   alloc(this->l.idx, x.lidx, this->l.size); alloc(this->l.start, x.lstart, this->l.startSize); alloc(this->l.row, x.lrow, this->l.startSize);
   alloc(this->l.ridx, x.lridx); alloc(this->l.rbeg, x.lrbeg);
   alloc(this->l.rorig, x.lrorig); alloc(this->l.rperm, x.lrperm);
   const auto copy = [](int* out, const std::vector<int>& in)
   {
      if(!in.empty()) std::memcpy(out, in.data(), in.size() * sizeof(int));
   };
   this->diag = x.diag; this->u.row.val = x.urval; this->u.col.val = x.ucval;
   this->l.val = x.lval; this->l.rval = x.lrval;
   copy(this->row.perm, x.rp); copy(this->row.orig, x.ro);
   copy(this->col.perm, x.cp); copy(this->col.orig, x.co);
   copy(this->u.row.idx, x.uridx); copy(this->u.row.start, x.urstart);
   copy(this->u.row.len, x.urlen); copy(this->u.row.max, x.urmax);
   copy(this->u.col.idx, x.ucidx); copy(this->u.col.start, x.ucstart);
   copy(this->u.col.len, x.uclen); copy(this->u.col.max, x.ucmax);
   copy(this->l.idx, x.lidx); copy(this->l.start, x.lstart); copy(this->l.row, x.lrow);
   copy(this->l.ridx, x.lridx); copy(this->l.rbeg, x.lrbeg);
   copy(this->l.rorig, x.lrorig); copy(this->l.rperm, x.lrperm);
   // Empty optional L rows must release old storage, as setupRowVals expects.
   if(x.lrbeg.empty())
   {
      if(this->l.ridx) spx_free(this->l.ridx);
      if(this->l.rbeg) spx_free(this->l.rbeg);
      if(this->l.rorig) spx_free(this->l.rorig);
      if(this->l.rperm) spx_free(this->l.rperm);
   }
   auto ring = [](typename CLUFactor<R>::Dring& list, typename CLUFactor<R>::Dring* elem,
                  const std::vector<int>& order)
   {
      auto* prev = &list;
      for(int idx : order)
      {
         auto* next = &elem[idx]; next->idx = idx;
         prev->next = next; next->prev = prev; prev = next;
      }
      prev->next = &list; list.prev = prev;
   };
   ring(this->u.row.list, this->u.row.elem, x.urorder);
   ring(this->u.col.list, this->u.col.elem, x.ucorder);
   this->u.row.size = int(x.uridx.size()); this->u.row.used = x.urused;
   this->u.col.size = int(x.ucidx.size()); this->u.col.used = x.ucused;
   this->l.size = int(x.lidx.size()); this->l.startSize = int(x.lstart.size());
   this->l.firstUpdate = x.firstUpdate; this->l.firstUnused = x.firstUnused;
   this->l.updateType = x.updateType; this->nzCnt = x.nz;
   this->initMaxabs = x.initMaxabs; this->maxabs = x.maxabs;
   this->clearLrightLookup(); this->stat = this->OK;
}

template <class R>
bool SLUFactor<R>::factorCacheRestore(int slot, const SVectorBase<R>* matrix[])
{
   if(slot < 0 || !factorSlots[slot].image) return false;
   bool started = false;
   try
   {
      std::vector<uint64_t> key;
      if(!factorCacheKey(matrix, key) || key != factorSlots[slot].key) return false;
      this->factorTime->start(); started = true;
      restoreFactorImage(*factorSlots[slot].image);
      this->factorTime->stop(); started = false;
      ++this->factorCount; ++factorHits;
      return true;
   }
   catch(const std::bad_alloc&) {}
   catch(const SPxMemoryException&) {}
   if(started) this->factorTime->stop();
   return false;
}

template <class R>
void SLUFactor<R>::factorCacheRemember(int slot, const SVectorBase<R>* matrix[])
{
   if(slot < 0 || factorSlots[slot].visits < 2) return;
   try
   {
      std::vector<uint64_t> key;
      if(!factorCacheKey(matrix, key)) return;
      auto image = factorImage();
      if(!image) return;
      int count = 0, oldest = -1;
      for(int i = 0; i < 4; ++i)
         if(i != slot && factorSlots[i].image)
         {
            ++count;
            if(oldest < 0 || factorSlots[i].tick < factorSlots[oldest].tick) oldest = i;
         }
      if(count >= 2) { factorSlots[oldest].image.reset(); factorSlots[oldest].key.clear(); }
      factorSlots[slot].key.swap(key);
      factorSlots[slot].image = std::move(image);
   }
   catch(const std::bad_alloc&) {}
}
} // namespace soplex
#endif
