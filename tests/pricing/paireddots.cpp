// SPDX-License-Identifier: Apache-2.0

#include "soplex.h"
#include "soplex/dsvector.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <random>
#include <stdexcept>

static void check(bool condition)
{
   if(!condition)
      throw std::runtime_error("paired pricing dots disagree with scalar StableSum dots");
}

static std::uint64_t bits(double value)
{
   std::uint64_t result;
   std::memcpy(&result, &value, sizeof(result));
   return result;
}

static double value(std::uint64_t raw)
{
   double result;
   std::memcpy(&result, &raw, sizeof(result));
   return result;
}

static void checkResult(double scalar, double paired)
{
   // IEEE NaN payload propagation is not portable across compilers. All other
   // results, including signed zeros and infinities, must match bit for bit.
   check(std::isnan(scalar) ? std::isnan(paired) : bits(scalar) == bits(paired));
}

int main()
{
   try
   {
      std::mt19937_64 random(223492347);
      int firstIndices[64], secondIndices[64];
      std::iota(firstIndices, firstIndices + 64, 0);
      std::iota(secondIndices, secondIndices + 64, 0);

      for(int trial = 0; trial < 4000; ++trial)
      {
         const int count = trial % 2 ? random() % 4 : random() % 65;
         std::shuffle(firstIndices, firstIndices + 64, random);
         std::shuffle(secondIndices, secondIndices + 64, random);
         soplex::DSVectorBase<double> first(64), second(64);
         soplex::VectorBase<double> work(64);
         auto number = [&]()
         {
            if(trial % 4 == 0)
               return double(int(random() % 2001) - 1000) / 8;
            std::uint64_t raw = random();
            if(trial % 4 == 2 && random() % 4 == 0)
               raw = random() & 1 ? 0 : UINT64_C(0x8000000000000000);
            // Mostly finite inputs, with nonfinite dot products possible by
            // overflow. Every fourth trial also permits input NaNs/infinities.
            if(trial % 4 != 3 && ((raw >> 52) & 2047) == 2047)
               raw ^= UINT64_C(1) << 52;
            return value(raw);
         };
         for(int k = 0; k < 64; ++k)
            work[k] = number();
         for(int k = 0; k < count; ++k)
         {
            // Keep explicit zeros to exercise the full element order.
            first.add(firstIndices[k]);
            first.value(k) = number();
            second.add(secondIndices[k]);
            second.value(k) = number();
         }
         double firstDot = 17, secondDot = 19;
         const bool handled = soplex::steeppr::pairedWeightDots(first, second, work,
                              firstDot, secondDot);
#if defined(__SSE2__) && !defined(SOPLEX_CHECK_STABLESUM)
         check(handled);
#else
         check(!handled);
#endif
         if(handled)
         {
            checkResult(first * work, firstDot);
            checkResult(second * work, secondDot);
         }
         if(count > 0)
         {
            second.remove(count - 1);
            firstDot = 17;
            secondDot = 19;
            check(!soplex::steeppr::pairedWeightDots(first, second, work, firstDot, secondDot));
            check(firstDot == 17 && secondDot == 19);
         }
      }

      // Generic scalar types retain their original dot operation.
      soplex::DSVectorBase<long double> column(1);
      soplex::VectorBase<long double> work(1);
      long double firstDot = 17, secondDot = 19;
      check(!soplex::steeppr::pairedWeightDots(column, column, work, firstDot, secondDot));
      check(firstDot == 17 && secondDot == 19);
      return 0;
   }
   catch(const std::exception& error)
   {
      std::fprintf(stderr, "%s\n", error.what());
      return 1;
   }
}
