// SPDX-License-Identifier: Apache-2.0
#include "soplex.h"
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

static void check(bool ok, const char* message)
{
   if(!ok) throw std::runtime_error(message);
}
static std::uint64_t bits(double value)
{
   std::uint64_t result;
   std::memcpy(&result, &value, sizeof(result));
   return result;
}
int main()
{
   try
   {
      std::vector<double> test(35, 1), weight(35, 2);
#if defined(__SSE2__) && defined(__FLT_EVAL_METHOD__) && __FLT_EVAL_METHOD__ == 0
      auto select = [&](double initial, double tolerance, int expectedIndex, double expected)
      {
         int index = -1;
         double best = initial;
         check(soplex::devexpr::densePriceCandidates(test.data(), weight.data(), 0, 35, 1,
               tolerance, best, index), "Real dense scan did not use candidate screening");
         check(index == expectedIndex, "Devex changed the first strictly best candidate");
         check(bits(best) == bits(expected), "Devex changed the winning price bits");
      };
      // The earlier equal maximum is in the second pair of a previous block.
      test[17] = test[18] = -4;
      select(0, 1e-6, 17, 8);
      select(8, 1e-6, -1, 8);
      test[17] = test[18] = 1;
      test[15] = test[16] = -6;
      weight[15] = weight[16] = 3;
      select(0, 1e-6, 15, 12);
      // The scalar remainder can improve the earlier block winner.
      test[34] = -5;
      weight[34] = .5;
      select(0, 1e-6, 34, 50);
      // An unordered incumbent never accepts a candidate or changes its payload.
      const std::uint64_t raw = UINT64_C(0x7ff8000000005a5a);
      double nan;
      std::memcpy(&nan, &raw, sizeof(nan));
      select(nan, 1e-6, -1, nan);
      test.assign(35, 1);
      weight.assign(35, 2);
      test[5] = test[6] = test[7] = test[8] = -1;
      weight[5] = nan;
      weight[6] = -0.0;
      weight[7] = weight[8] = 0.0;
      select(0, 0, 7, std::numeric_limits<double>::infinity());
      select(std::numeric_limits<double>::infinity(), 0, -1,
             std::numeric_limits<double>::infinity());
      // A weight below the floor uses the same denominator as scalar pricing.
      test.assign(35, 1);
      weight.assign(35, 2);
      test[9] = -2;
      weight[9] = -8;
      select(0, .5, 9, 8);
#else
      int unused = -1;
      double incumbent = 0;
      check(!soplex::devexpr::densePriceCandidates(test.data(), weight.data(), 0, 35, 1,
            1e-6, incumbent, unused), "Unsupported evaluation did not retain scalar pricing");
#endif
      double best = 1;
      int index = 17;
      check(!soplex::devexpr::densePriceCandidates(test.data(), weight.data(), 0, 15, 1,
            1e-6, best, index), "Short scan did not retain scalar pricing");
      check(!soplex::devexpr::densePriceCandidates(test.data(), weight.data(), 0, 35, 2,
            1e-6, best, index), "Strided scan did not retain scalar pricing");
      check(best == 1 && index == 17, "Declined scan changed the incumbent");
      long double generic[16] = {}, genericBest = 1;
      check(!soplex::devexpr::densePriceCandidates(generic, generic, 0, 16, 1,
            1e-6L, genericBest, index), "Generic type did not retain scalar pricing");
      check(genericBest == 1 && index == 17, "Generic fallback changed the incumbent");
      return 0;
   }
   catch(const std::exception& error)
   {
      std::fprintf(stderr, "%s\n", error.what());
      return 1;
   }
}
