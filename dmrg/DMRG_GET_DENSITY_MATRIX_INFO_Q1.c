#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "dmrg.h"

DMRG_DENSITY_MATRIX_INFO *DMRG_GET_DENSITY_MATRIX_INFO_Q1(int *Q_Number_LLLR, int dim_LLLR) {
   
   //Check Point
   if (dim_LLLR <= 0) {
      printf("Error in DMRG_GET_DENSITY_MATRIX_INFO_Q1\n");
      printf("dim_LLLR = %d\n", dim_LLLR);
      exit(1);
   }
   
   int basis;
   int block_num = 1;
   
   for (basis = 0; basis < dim_LLLR - 1; basis++) {
      //Check Point
      if (Q_Number_LLLR[basis] > Q_Number_LLLR[basis + 1]) {
         printf("Error in DMRG_GET_DENSITY_MATRIX_INFO_Q1\n");
         printf("Q_Number_LLLR[%d] (%d) > Q_Number_LLLR[%d] (%d)\n",basis, Q_Number_LLLR[basis], basis + 1, Q_Number_LLLR[basis + 1]);
         exit(1);
      }
      if (Q_Number_LLLR[basis] != Q_Number_LLLR[basis + 1]) {
         block_num++;
      }
   }
   
   //Allocate
   DMRG_DENSITY_MATRIX_INFO *Info = malloc(sizeof(DMRG_DENSITY_MATRIX_INFO));
   Info->Dim_Block       = GET_ARRAY_INT1(block_num);
   Info->Q_Number1_Block = GET_ARRAY_INT1(block_num);
   Info->Q_Number2_Block = GET_ARRAY_INT1(1);
   Info->Q_Number3_Block = GET_ARRAY_INT1(1);
   Info->block_num       = block_num;
   
   int min_q_num = Q_Number_LLLR[0];
   int max_q_num = Q_Number_LLLR[dim_LLLR - 1];
   int *Temp_Dim_Block = GET_ARRAY_INT1(max_q_num - min_q_num + 1);
   
   for (basis = 0; basis < dim_LLLR; basis++) {
      Temp_Dim_Block[Q_Number_LLLR[basis] - min_q_num]++;
   }
   
   int q_num;
   block_num = 0;
   for (q_num = 0; q_num <= max_q_num - min_q_num; q_num++) {
      if (Temp_Dim_Block[q_num] > 0) {
         Info->Dim_Block[block_num]       = Temp_Dim_Block[q_num];
         Info->Q_Number1_Block[block_num] = q_num + min_q_num;
         block_num++;
      }
   }
   
   //Check Point
   if (block_num != Info->block_num) {
      printf("Error in DMRG_GET_DENSITY_MATRIX_INFO_Q1\n");
      printf("%d != %d\n", block_num, Info->block_num);
      exit(1);
   }
   
   FREE_ARRAY_INT1(Temp_Dim_Block);
   
   return Info;
   
}
