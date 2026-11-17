// SPDX-License-Identifier: Apache-2.0

#include "soplex/slufactor.h"
#include "soplex/dsvector.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <vector>

static void check(bool condition)
{
   if(!condition)
      throw std::runtime_error("Sparse right solve disagrees with the known solution");
}

int main()
{
   try
   {
      const int maxdim = 16384;
      auto tolerances = std::make_shared<soplex::Tolerances>();
      soplex::SLUFactor<double> factor;
      factor.setTolerances(tolerances);
      std::vector<soplex::DSVectorBase<double>> matrix;
      matrix.reserve(maxdim);
      std::vector<const soplex::SVectorBase<double>*> columns(maxdim);
      for(int c = 0; c < maxdim; ++c)
      {
         matrix.emplace_back(2);
         columns[c] = &matrix[c];
      }
      // Change row permutation and dimension, then grow the basis again.
      // The marker buffer's capacity and lookup validity must follow each load.
      for(int reload = 0; reload < 3; ++reload)
      {
         const int dim = reload == 1 ? maxdim / 2 : maxdim;
         for(int c = 0; c < dim; ++c)
         {
            matrix[c].clear();
            const int first = (c / 2 * 2 + reload * 13) % dim;
            matrix[c].add(first, c % 2 ? 1.0 : 8.0);
            matrix[c].add((first + 1) % dim, c % 2 ? 8.0 : 1.0);
         }
         check(factor.load(columns.data(), dim) == soplex::SLinSolver<double>::OK);
         for(int trial = 0; trial < 24; ++trial)
         {
            std::vector<double> expected(dim, 0.0), expected2(dim, 0.0);
            std::vector<double> denseRhs(dim, 0.0), denseRhs2(dim, 0.0);
            for(int j = 0; j < 8; ++j)
            {
               expected[((j * 137 + trial * 47) % (dim / 2)) * 2] = 0.25 * (j + 1);
               expected2[((j * 61 + trial * 127) % (dim / 2)) * 2 + 1] = -0.125 * (j + 1);
            }
            for(int c = 0; c < dim; ++c)
               for(int k = 0; k < matrix[c].size(); ++k)
               {
                  const int row = matrix[c].index(k);
                  denseRhs[row] += expected[c] * matrix[c].value(k);
                  denseRhs2[row] += expected2[c] * matrix[c].value(k);
               }
            soplex::DSVectorBase<double> rhs(dim), rhs2(dim);
            for(int row = 0; row < dim; ++row)
            {
               if(denseRhs[row] != 0)
                  rhs.add(row, denseRhs[row]);
               if(denseRhs2[row] != 0)
                  rhs2.add(row, denseRhs2[row]);
            }
            soplex::SSVectorBase<double> single(dim, tolerances), x(dim, tolerances);
            soplex::SSVectorBase<double> y(dim, tolerances), secondRhs(dim, tolerances);
            factor.solveRight4update(single, rhs);
            secondRhs = rhs2;
            factor.solve2right4update(x, y, rhs, secondRhs);
            for(int row = 0; row < dim; ++row)
            {
               check(std::isfinite(single[row]) && std::abs(single[row] - expected[row]) < 1e-10);
               check(std::isfinite(x[row]) && std::abs(x[row] - expected[row]) < 1e-10);
               check(std::isfinite(y[row]) && std::abs(y[row] - expected2[row]) < 1e-10);
            }
         }
      }
      return 0;
   }
   catch(const std::exception& error)
   {
      std::fprintf(stderr, "%s\n", error.what());
      return 1;
   }
}
