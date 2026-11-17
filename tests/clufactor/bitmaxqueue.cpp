// SPDX-License-Identifier: Apache-2.0

#include "soplex/clufactor.h"

#include <algorithm>
#include <cstdio>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <vector>

static void check(bool condition)
{
   if(!condition)
      throw std::runtime_error("LU integer queue disagrees with reference maximum queue");
}

int main()
{
   try
   {
      std::mt19937 random(1703);
      const int dimensions[] = {1, 2, 63, 64, 65, 127, 128, 129, 4095, 4096, 4097, 65535, 65536};
      for(int dim : dimensions)
      {
         std::vector<int> indices(dim);
         std::iota(indices.begin(), indices.end(), 0);
         std::shuffle(indices.begin(), indices.end(), random);
         std::priority_queue<int> reference;
         std::vector<unsigned char> pending(dim, 0);
         int count = dim / 2;
         for(int i = 0; i < count; ++i)
         {
            reference.push(indices[i]);
            pending[indices[i]] = 1;
         }
         soplex::CLUFactorBitMaxQueue queue(indices.data(), count);
         for(int step = 0; step < 50000; ++step)
         {
            if(reference.empty() || (count < dim / 2 && random() % 3 != 0))
            {
               int key = random() % dim;
               while(pending[key])
                  key = (key + 1) % dim;
               queue.push(key, &count);
               reference.push(key);
               pending[key] = 1;
            }
            else
            {
               check(queue.maximum() == reference.top());
               int key = queue.pop(&count);
               check(key == reference.top());
               reference.pop();
               pending[key] = 0;
            }
            check(count == int(reference.size()));
         }
         while(!reference.empty())
         {
            check(queue.maximum() == reference.top());
            check(queue.pop(&count) == reference.top());
            reference.pop();
         }
         check(count == 0);
         // Reuse inactive groups and words after the complete queue has emptied.
         for(int key = 0; key < dim; ++key)
            queue.push(key, &count);
         for(int key = dim - 1; key >= 0; --key)
         {
            check(queue.maximum() == key);
            check(queue.pop(&count) == key);
         }
         check(count == 0);
      }
   }
   catch(const std::exception& error)
   {
      std::fprintf(stderr, "%s\n", error.what());
      return 1;
   }
   std::puts("LU integer queue: all word/summary boundaries, mixed operations, and reuse pass");
   return 0;
}
