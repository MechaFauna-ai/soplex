/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxdefines.hpp
 * @brief General templated functions for SoPlex
 */

// Defining the static members of the Param class

#include <cstdint>
#include <cstring>
#include <type_traits>

namespace soplex
{

/// returns \c true iff |a-b| <= eps
template <class R, class S, class T>
inline bool EQ(R a, S b, T eps)
{
   return spxAbs(R(a - b)) <= eps;
}

/// returns \c true iff |a-b| > eps
template <class R, class S, class T>
inline bool NE(R a, S b, T eps)
{
   return spxAbs(a - b) > eps;
}

/// returns \c true iff b-a > eps
template <class R, class S, class T>
inline bool LT(R a, S b, T eps)
{
   return (a - b) < -eps;
}

/// returns \c true iff a-b <= eps
template <class R, class S, class T>
inline bool LE(R a, S b, T eps)
{
   return (a - b) <= eps;
}

/// returns \c true iff a-b > eps
template <class R, class S, class T>
inline bool GT(R a, S b, T eps)
{
   return (a - b) > eps;
}

/// returns \c true iff b-a <= eps
template <class R, class S, class T>
inline bool GE(R a, S b, T eps)
{
   return (a - b) >= -eps;
}

/// returns \c true iff |a| <= eps
template <class R, class T>
inline bool isZero(R a, T eps)
{
   return spxAbs(a) <= eps;
}

/// returns \c true iff |a| > eps
template <class R, class T>
inline bool isNotZero(R a, T eps)
{
   return spxAbs(a) > eps;
}

/// returns \c true iff |relDiff(a,b)| <= eps
template <class R, class S, class T>
inline bool EQrel(R a, S b, T eps)
{
   return spxAbs(relDiff(a, b)) <= eps;
}

/// returns \c true iff |relDiff(a,b)| > eps
template <class R, class S, class T>
inline bool NErel(R a, S b, T eps)
{
   return spxAbs(relDiff(a, b)) > eps;
}

/// returns \c true iff relDiff(b,a) > eps
template <class R, class S, class T>
inline bool LTrel(R a, S b, T eps)
{
   return relDiff(a, b) < -eps;
}

/// returns \c true iff relDiff(a,b) <= eps
template <class R, class S, class T>
inline bool LErel(R a, S b, T eps)
{
   return relDiff(a, b) <= eps;
}

/// returns \c true iff relDiff(a,b) > eps
template <class R, class S, class T>
inline bool GTrel(R a, S b, T eps)
{
   return relDiff(a, b) > eps;
}

/// returns \c true iff relDiff(b,a) <= eps
template <class R, class S, class T>
inline bool GErel(R a, S b, T eps)
{
   return relDiff(a, b) >= -eps;
}

// Instantiation for Real
inline Real spxLdexp(Real x, int exp)
{
   // fast path for double: if x is a normal number and x * 2^exp is a normal number, too, the result is obtained exactly
   // by adding exp to the exponent field, which is what ldexp() computes; all other cases (zero, subnormal, infinite,
   // NaN, overflow, underflow) are handled by ldexp()
   if(std::is_same<Real, double>::value && exp > -2048 && exp < 2048)
   {
      const double d = double(x);
      uint64_t bits;
      std::memcpy(&bits, &d, sizeof(uint64_t));
      const int biasedexp = int((bits >> 52) & 0x7ff);
      const int newexp = biasedexp + exp;

      if(biasedexp != 0 && biasedexp != 0x7ff && newexp > 0 && newexp < 0x7ff)
      {
         double result;
         bits = (bits & ~(uint64_t(0x7ff) << 52)) | (uint64_t(newexp) << 52);
         std::memcpy(&result, &bits, sizeof(uint64_t));
         assert(result == std::ldexp(d, exp));
         return Real(result);
      }
   }

   return std::ldexp(x, exp);
}

inline Real spxFrexp(Real x, int* exp)
{
   return std::frexp(x, exp);
}

} // namespace soplex
