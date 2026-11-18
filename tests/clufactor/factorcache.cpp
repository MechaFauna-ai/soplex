// SPDX-License-Identifier: Apache-2.0
#include "soplex/slufactor.h"
#include "soplex/dsvector.h"
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <random>
#include <stdexcept>
#include <vector>
using namespace soplex;
static void check(bool ok)
{
   if(!ok) throw std::runtime_error("Cached LU differs from a fresh factorization");
}
class Probe : public SLUFactor<double>
{
public:
   std::vector<uint64_t> words() const
   {
      std::vector<uint64_t> out;
      auto integer = [&](int v) { out.push_back(uint64_t(v)); };
      auto value = [&](double v) { uint64_t bits; std::memcpy(&bits, &v, sizeof(bits)); out.push_back(bits); };
      integer(this->stat); integer(this->nzCnt); value(this->initMaxabs); value(this->maxabs);
      for(int i = 0; i < this->thedim; ++i)
      {
         integer(this->row.perm[i]); integer(this->row.orig[i]);
         integer(this->col.perm[i]); integer(this->col.orig[i]); value(this->diag[i]);
         integer(this->u.row.len[i]); integer(this->u.col.len[i]);
         for(int j = this->u.row.start[i]; j < this->u.row.start[i] + this->u.row.len[i]; ++j)
         { integer(this->u.row.idx[j]); value(this->u.row.val[j]); }
         for(int j = this->u.col.start[i]; j < this->u.col.start[i] + this->u.col.len[i]; ++j)
         { integer(this->u.col.idx[j]); value(this->u.col.val[j]); }
#ifdef SOPLEX_WITH_L_ROWS
         integer(this->l.rorig[i]); integer(this->l.rperm[i]);
         integer(this->l.rbeg[i + 1] - this->l.rbeg[i]);
         for(int j = this->l.rbeg[i]; j < this->l.rbeg[i + 1]; ++j)
         { integer(this->l.ridx[j]); value(this->l.rval[j]); }
#endif
      }
      integer(this->l.firstUpdate); integer(this->l.firstUnused);
      for(int i = 0; i < this->l.firstUnused; ++i)
      { integer(this->l.row[i]); integer(this->l.start[i + 1] - this->l.start[i]); }
      for(int i = 0; i < this->l.start[this->l.firstUnused]; ++i)
      { integer(this->l.idx[i]); value(this->l.val[i]); }
      return out;
   }
};
static void same(const VectorBase<double>& a, const VectorBase<double>& b)
{
   check(a.dim() == b.dim());
   for(int i = 0; i < a.dim(); ++i) check(std::memcmp(&a[i], &b[i], sizeof(double)) == 0);
}
int main()
{
   std::mt19937 rng(5813);
   uint64_t hits = 0;
   for(int mode = 0; mode < 2; ++mode)
   {
      Probe cached;
      auto tol = std::make_shared<Tolerances>();
      cached.setTolerances(tol);
      auto updateType = mode ? SLUFactor<double>::ETA : SLUFactor<double>::FOREST_TOMLIN;
      cached.setUtype(updateType);
      for(int dim : {128, 160, 96, 128})
      {
         std::vector<DSVectorBase<double>> a(dim), b(dim);
         std::vector<const SVectorBase<double>*> A(dim), B(dim);
         for(int c = 0; c < dim; ++c)
         {
            for(int r = 0; r < dim; ++r)
               if(r == c || rng() % 100 < 20)
               {
                  double v = r == c ? dim : (int(rng() % 17) - 8) / 16.0;
                  if(v != 0) { a[c].add(r, v); b[c].add(r, v + (r == c ? 1.0 : 0.0)); }
               }
            A[c] = &a[c]; B[c] = &b[c];
         }
         VectorBase<double> rhs(dim), x(dim), y(dim);
         for(int r = 0; r < dim; ++r) rhs[r] = int(rng() % 17) - 8;
         for(int step = 0; step < 24; ++step)
         {
            if(step % 3 == 0) cached.clear();
            if(step == 6) a[0].value(0) += 0.25;
            if(step == 12) a[0].value(0) -= 0.25;
            if(step == 8) tol->setEpsilonFactorization(2e-16);
            if(step == 14) tol->setEpsilonFactorization(SOPLEX_DEFAULT_EPS_FACTOR);
            if(step == 10) cached.setMarkowitz(0.1);
            if(step == 16) cached.setMarkowitz(0.01);
            auto columns = step % 2 ? B.data() : A.data();
            check(cached.load(columns, dim) == SLinSolver<double>::OK);
            // A new solver cannot have an admitted image: the reference always factors.
            Probe fresh;
            fresh.setTolerances(tol); fresh.setUtype(updateType);
            fresh.setMarkowitz(cached.markowitz());
            check(fresh.load(columns, dim) == SLinSolver<double>::OK);
            check(cached.words() == fresh.words());
            cached.solveRight(x, rhs); fresh.solveRight(y, rhs); same(x, y);
            cached.solveLeft(x, rhs); fresh.solveLeft(y, rhs); same(x, y);
            DSVectorBase<double> replacement(*columns[0]);
            replacement.value(0) += 0.125;
            SSVectorBase<double> update(dim, tol), updateRef(dim, tol);
            cached.solveRight4update(update, replacement);
            fresh.solveRight4update(updateRef, replacement);
            for(int i = 0; i < dim; ++i) check(update[i] == updateRef[i]);
            check(cached.change(0, replacement) == SLinSolver<double>::OK);
            check(fresh.change(0, replacement) == SLinSolver<double>::OK);
            check(cached.words() == fresh.words());
            cached.solveRight(x, rhs); fresh.solveRight(y, rhs); same(x, y);
            cached.solveLeft(x, rhs); fresh.solveLeft(y, rhs); same(x, y);
         }
         const auto before = cached.factorCacheHits();
         check(std::fesetround(FE_DOWNWARD) == 0);
         check(cached.load(A.data(), dim) == SLinSolver<double>::OK);
         check(std::fesetround(FE_TONEAREST) == 0);
         check(cached.factorCacheHits() == before);
      }
      hits += cached.factorCacheHits();
   }
#if defined(__SSE2__) && FLT_EVAL_METHOD == 0
   if(std::is_same<Real, double>::value) check(hits > 0);
#endif
}
