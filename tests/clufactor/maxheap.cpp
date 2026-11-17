// SPDX-License-Identifier: Apache-2.0

#include "soplex/clufactor.h"

#include <algorithm>
#include <cstdio>
#include <functional>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <vector>

static void check(bool condition)
{
   if(!condition)
      throw std::runtime_error("LU heap disagrees with its priority queue");
}

int main()
{
   try
   {
      std::mt19937 random(4271);
      for(int trial = 0; trial < 200; ++trial)
      {
         const int dim = 1 + random() % 4096;
         std::vector<int> heap(dim), permutation(dim);
         std::iota(heap.begin(), heap.end(), 0);
         std::shuffle(heap.begin(), heap.end(), random);
         for(int i = 0; i < dim; ++i)
            permutation[i] = trial % 2 ? int(random() % 64) : i;
         std::priority_queue<int> reference;
         for(int i : heap)
            reference.push(permutation[i]);
         soplex::initMaxHeap(heap.data(), dim, permutation.data());
         for(int i = 1; i < dim; ++i)
            check(heap[(i - 1) / 4] >= heap[i]);
         int count = dim;
         for(int step = 0; step < 2 * dim; ++step)
         {
            if(count == 0 || (count < dim && random() % 3 != 0))
            {
               const int key = int(random() % 127);
               soplex::enQueueMax(heap.data(), &count, key);
               reference.push(key);
            }
            else
            {
               check(soplex::deQueueMax(heap.data(), &count) == reference.top());
               reference.pop();
            }
            check(count == int(reference.size()));
         }
         while(count > 0)
         {
            check(soplex::deQueueMax(heap.data(), &count) == reference.top());
            reference.pop();
         }

         // Generic construction retains binary layout for min-heap consumers,
         // including transparent and lambda comparator types.
         for(int kind = 0; kind < 3; ++kind)
         {
            std::iota(heap.begin(), heap.end(), 0);
            if(kind == 0)
               soplex::initHeap(heap.data(), dim, permutation.data(), std::greater<int>());
            else if(kind == 1)
               soplex::initHeap(heap.data(), dim, permutation.data(), std::greater<>());
            else
               soplex::initHeap(heap.data(), dim, permutation.data(), [](int a, int b) { return a > b; });
            for(int i = 1; i < dim; ++i)
               check(heap[(i - 1) / 2] <= heap[i]);
            std::priority_queue<int, std::vector<int>, std::greater<int>> minimum;
            for(int value : permutation)
               minimum.push(value);
            count = dim;
            while(count > 0)
            {
               check(soplex::deQueueMin(heap.data(), &count) == minimum.top());
               minimum.pop();
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
