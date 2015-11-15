#ifndef LEVI_CIVITA_H
#define LEVI_CIVITA_H

#include "tensor.h"
#include <vector>
#include <algorithm>


namespace radmat
{


  template<idx_t RANK>
    int bubble_sort(const idx_t input_sequence[])
    {

      // bubble sort will destroy this guy doing the permutations
      idx_t sequence[RANK];
      std::copy(input_sequence,input_sequence+RANK,sequence); 


      idx_t i, j, temp, nperms;

      nperms = 0; 

      for (i = (RANK - 1); i > 0; i--)
      {
        for (j = 1; j <= i; j++)
        {
          if (sequence[j-1] > sequence[j])
          {
            ++nperms;
            temp = sequence[j-1];
            sequence[j-1] = sequence[j];
            sequence[j] = temp;
          }
        }
      }

      return nperms; 
    }


  template<typename T, idx_t RANK>
    Tensor<T,RANK> levi_civita(void)  // eps^{0...RANK} = 1
    {

      std::vector<idx_t> shape(RANK,RANK);
      Tensor<T,RANK> levi(shape, T(0)); 
      idx_t sequence[RANK];

      for(idx_t i = 0; i < RANK; ++i)
        sequence[i] = i;

      T plus(1);
      T minus(-1);

      do
      {

        int nperm = bubble_sort<RANK>(sequence);
        levi.getElem(sequence) = (nperm % 2 == 0) ? plus : minus; 

      } while (std::next_permutation(sequence,sequence+RANK));


      return levi; 
    }



}

#endif /* LEVI_CIVITA_H */
