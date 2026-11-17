// SPDX-License-Identifier: Apache-2.0
#include "soplex/slufactor.h"
#include "soplex/dsvector.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>
static void check(bool ok)
{
   if(!ok)
      throw std::runtime_error("LU column growth changed the known solution");
}
int main()
{
   try
   {
      std::mt19937 rng(4173);
      auto tolerances = std::make_shared<soplex::Tolerances>();
      for(int mode = 0; mode < 2; ++mode)
      {
         soplex::SLUFactor<double> factor;
         factor.setTolerances(tolerances);
         factor.setUtype(mode ? soplex::SLUFactor<double>::ETA : soplex::SLUFactor<double>::FOREST_TOMLIN);
         for(int dim : {128, 96, 192})
         {
            std::vector<soplex::DSVectorBase<double>> matrix;
            matrix.reserve(dim);
            std::vector<const soplex::SVectorBase<double>*> columns(dim);
            for(int c = 0; c < dim; ++c)
            {
               matrix.emplace_back(dim);
               columns[c] = &matrix[c];
               for(int row = 0; row < dim; ++row)
                  if(row == c || rng() % 100 < 12)
                  {
                     const double value = row == c ? dim : (int(rng() % 17) - 8) / 8.0;
                     if(value != 0)
                        matrix[c].add(row, value);
                  }
            }
            check(factor.load(columns.data(), dim) == soplex::SLinSolver<double>::OK);
            const int changes = mode ? 1 : 8;
            for(int step = 0; step <= changes; ++step)
            {
               std::vector<double> expected(dim), right(dim, 0.0), left(dim, 0.0);
               for(int row = 0; row < dim; ++row)
                  expected[row] = (int(rng() % 17) - 8) / 8.0;
               for(int c = 0; c < dim; ++c)
                  for(int k = 0; k < matrix[c].size(); ++k)
                  {
                     right[matrix[c].index(k)] += matrix[c].value(k) * expected[c];
                     left[c] += matrix[c].value(k) * expected[matrix[c].index(k)];
                  }
               soplex::VectorBase<double> rhs(dim), rhsLeft(dim), result(dim);
               soplex::DSVectorBase<double> sparseRhs(dim);
               for(int row = 0; row < dim; ++row)
               {
                  rhs[row] = right[row];
                  rhsLeft[row] = left[row];
                  if(right[row] != 0)
                     sparseRhs.add(row, right[row]);
               }
               factor.solveRight(result, rhs);
               for(int row = 0; row < dim; ++row)
                  check(std::isfinite(result[row]) && std::abs(result[row] - expected[row]) < 1e-10);
               factor.solveLeft(result, rhsLeft);
               for(int row = 0; row < dim; ++row)
                  check(std::isfinite(result[row]) && std::abs(result[row] - expected[row]) < 1e-10);
               soplex::SSVectorBase<double> sparseResult(dim, tolerances);
               factor.solveRight4update(sparseResult, sparseRhs);
               for(int row = 0; row < dim; ++row)
                  check(std::isfinite(sparseResult[row]) && std::abs(sparseResult[row] - expected[row]) < 1e-10);
               if(step == changes)
                  break;
               const int c = step * 31 % dim;
               soplex::DSVectorBase<double> replacement(dim);
               for(int row = 0; row < dim; ++row)
                  if(row == c || rng() % 100 < 70)
                  {
                     const double value = row == c ? dim + step : (int(rng() % 17) - 8) / 8.0;
                     if(value != 0)
                        replacement.add(row, value);
                  }
               soplex::SSVectorBase<double> update(dim, tolerances);
               factor.solveRight4update(update, replacement);
               check(factor.change(c, replacement) == soplex::SLinSolver<double>::OK);
               matrix[c] = replacement;
            }
         }
      }
      std::puts("PASS LU column growth through factorization, FT/eta updates and dimension reloads");
      return 0;
   }
   catch(const std::exception& e)
   {
      std::fprintf(stderr, "%s\n", e.what());
      return 1;
   }
}
