#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <omp.h>
#include "dmrg.h"

DMRG_SYSTEM_INFO *DMRG_GET_SYSTEM_INFO_Q1(int *Q_Number, DMRG_BASIS *Dmrg_Basis, double *Vec, int max_dim, int p_threads, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->make_dens_mat = omp_get_wtime();
   
   int dim_LLLR = Dmrg_Basis->dim_LLLR;
   int block,iter1,iter2,count,dim;
   double val1,val2;
   
   DMRG_DENSITY_MATRIX_INFO *Info = DMRG_GET_DENSITY_MATRIX_INFO_Q1(Q_Number, dim_LLLR);
   
   double ***Density_Matrix = malloc(sizeof(double**)*Info->block_num);
   double ***Vector_DM      = malloc(sizeof(double**)*Info->block_num);
   double **Value_DM        = malloc(sizeof(double*)*Info->block_num);
   
   for (block = 0; block < Info->block_num; block++) {
      Density_Matrix[block] = DMRG_GET_DENSITY_MATRIX(block, Info->Dim_Block, Dmrg_Basis, Vec, p_threads);
      Vector_DM[block]      = GET_ARRAY_DOUBLE2(Info->Dim_Block[block], Info->Dim_Block[block]);
      Value_DM[block]       = GET_ARRAY_DOUBLE1(Info->Dim_Block[block]);
   }
   
#pragma omp parallel for private (dim) num_threads (p_threads)
   for (block = 0; block < Info->block_num; block++) {
      dim = Info->Dim_Block[block];
      LAPACK_DSYEV(Density_Matrix[block], dim, dim, Value_DM[block], Vector_DM[block], dim, dim);
   }
   
   int *Block_Sorted    = GET_ARRAY_INT1(dim_LLLR);
   int *Basis_Sorted    = GET_ARRAY_INT1(dim_LLLR);
   int *Q_Number_Sorted = GET_ARRAY_INT1(dim_LLLR);
   
   count = 0;
   for (block = 0; block < Info->block_num; block++) {
      for (iter1 = 0; iter1 < Info->Dim_Block[block]; iter1++) {
         Block_Sorted[count]    = block;
         Basis_Sorted[count]    = iter1;
         Q_Number_Sorted[count] = Info->Q_Number1_Block[block];
         count++;
      }
   }
   
   //Check Point 1
   if (count != dim_LLLR) {
      printf("Error in DMRG_MAKE_TRANSFORMATION_MATRIX_Q1 at Check Point 1\n");
      exit(1);
   }
   
   //Bubble Sort
   for (iter1 = 0; iter1 < dim_LLLR - 1; iter1++) {
      count = 0;
      for (iter2 = 0; iter2 < dim_LLLR - 1 - iter1; iter2++) {
         val1 = Value_DM[Block_Sorted[iter2]][Basis_Sorted[iter2]];
         val2 = Value_DM[Block_Sorted[iter2 + 1]][Basis_Sorted[iter2 + 1]];
         if (val1 < val2) {
            SWAP_INT(&Block_Sorted[iter2]   , &Block_Sorted[iter2 + 1]);
            SWAP_INT(&Basis_Sorted[iter2]   , &Basis_Sorted[iter2 + 1]);
            SWAP_INT(&Q_Number_Sorted[iter2], &Q_Number_Sorted[iter2 + 1]);
            count = 1;
         }
      }
      if (count == 0) {
         break;
      }
   }
   
   //Allocate
   DMRG_SYSTEM_INFO *Dmrg_System = malloc(sizeof(DMRG_SYSTEM_INFO));
   Dmrg_System->Val_DM_Dist = GET_ARRAY_DOUBLE1(dim_LLLR);
   Dmrg_System->sum_val_dm = 0;
   
   //Store System Information
   for (iter1 = 0; iter1 < dim_LLLR; iter1++) {
      Dmrg_System->Val_DM_Dist[iter1] = Value_DM[Block_Sorted[iter1]][Basis_Sorted[iter1]];
      Dmrg_System->sum_val_dm = Dmrg_System->sum_val_dm + Value_DM[Block_Sorted[iter1]][Basis_Sorted[iter1]];
   }
   
   int dim_renorm = 0;
   double zero = pow(10,-15);
   
   
   if (dim_LLLR <= max_dim) {
      for (iter1 = dim_LLLR; iter1 >= 1; iter1--) {
         val1 = Value_DM[Block_Sorted[iter1 - 1]][Basis_Sorted[iter1 - 1]];
         if (val1 > zero) {
            dim_renorm = iter1;
            break;
         }
      }
   }
   else {
      for (iter1 = max_dim; iter1 >= 1; iter1--) {
         val1 = Value_DM[Block_Sorted[iter1 - 1]][Basis_Sorted[iter1 - 1]];
         val2 = Value_DM[Block_Sorted[iter1]][Basis_Sorted[iter1]];
         if (fabs(val1 - val2) > zero && val1 > zero) {
            dim_renorm = iter1;
            break;
         }
      }
   }
   
   //Check Point 2
   if (dim_renorm <= 0) {
      printf("Error in DMRG_MAKE_TRANSFORMATION_MATRIX_Q1 at Check Point 2\n");
      exit(1);
   }
   
   //Store System Information: Truncation Error
   Dmrg_System->tr_error = 0;
   for (iter1 = 0; iter1 < dim_renorm; iter1++) {
      Dmrg_System->tr_error = Dmrg_System->tr_error + Value_DM[Block_Sorted[iter1]][Basis_Sorted[iter1]];
   }
   Dmrg_System->tr_error = Dmrg_System->sum_val_dm - Dmrg_System->tr_error;
   
   for (iter1 = 0; iter1 < dim_renorm - 1; iter1++) {
      count = 0;
      for (iter2 = 0; iter2 < dim_renorm - 1 - iter1; iter2++) {
         if (Q_Number_Sorted[iter2] > Q_Number_Sorted[iter2 + 1]) {
            SWAP_INT(&Block_Sorted[iter2]   , &Block_Sorted[iter2 + 1]);
            SWAP_INT(&Basis_Sorted[iter2]   , &Basis_Sorted[iter2 + 1]);
            SWAP_INT(&Q_Number_Sorted[iter2], &Q_Number_Sorted[iter2 + 1]);
            count = 1;
         }
      }
      if (count == 0) {
         break;
      }
   }
   
   int *Row_Base = GET_ARRAY_INT1(Info->block_num + 1);
   
   for (iter1 = 0; iter1 < Info->block_num; iter1++) {
      Row_Base[iter1 + 1] = Row_Base[iter1] + Info->Dim_Block[iter1];
   }
   
   long elem_num = 0;
   for (iter1 = 0; iter1 < dim_renorm; iter1++) {
      elem_num = elem_num + Info->Dim_Block[Block_Sorted[iter1]];
   }
   
   //Allocate
   Dmrg_System->Trans_Matrix        = GET_CCS1(dim_renorm, elem_num);
   Dmrg_System->Trans_Matrix_Dagger = GET_CRS1(dim_renorm, elem_num);
   Dmrg_System->Q_Number1           = GET_ARRAY_INT1(dim_renorm);
   Dmrg_System->Q_Number2           = GET_ARRAY_INT1(1);
   Dmrg_System->Q_Number3           = GET_ARRAY_INT1(1);
   Dmrg_System->dim_renorm          = dim_renorm;
   Dmrg_System->dim_LLLR            = dim_LLLR;
   
   //Store Quantum Number
   for (iter1 = 0; iter1 < dim_renorm; iter1++) {
      Dmrg_System->Q_Number1[iter1] = Q_Number_Sorted[iter1];
   }
   
   //Store Dummy Value
   Dmrg_System->Q_Number2[0] = INT_MIN;
   Dmrg_System->Q_Number3[0] = INT_MIN;
   
   elem_num = 0;
   for (iter1 = 0; iter1 < dim_renorm; iter1++) {
      for (iter2 = 0; iter2 < Info->Dim_Block[Block_Sorted[iter1]]; iter2++) {
         Dmrg_System->Trans_Matrix->Val[elem_num]        = Vector_DM[Block_Sorted[iter1]][Basis_Sorted[iter1]][iter2];
         Dmrg_System->Trans_Matrix->Row[elem_num]        = Row_Base[Block_Sorted[iter1]] + iter2;
         Dmrg_System->Trans_Matrix_Dagger->Val[elem_num] = Vector_DM[Block_Sorted[iter1]][Basis_Sorted[iter1]][iter2];
         Dmrg_System->Trans_Matrix_Dagger->Col[elem_num] = Row_Base[Block_Sorted[iter1]] + iter2;
         elem_num++;
      }
      Dmrg_System->Trans_Matrix->Col[iter1 + 1]        = elem_num;
      Dmrg_System->Trans_Matrix_Dagger->Row[iter1 + 1] = elem_num;
   }
   
   //Check Point 3
   if (elem_num != Dmrg_System->Trans_Matrix->max_val) {
      printf("Error in DMRG_MAKE_TRANSFORMATION_MATRIX_Q1 at Check Point 3\n");
      exit(1);
   }
   
   Dmrg_System->Trans_Matrix->row_dim        = dim_LLLR;
   Dmrg_System->Trans_Matrix->col_dim        = dim_renorm;
   Dmrg_System->Trans_Matrix_Dagger->row_dim = dim_renorm;
   Dmrg_System->Trans_Matrix_Dagger->col_dim = dim_LLLR;
   
   //FREE_ARRAY
   for (block = 0; block < Info->block_num; block++) {
      FREE_ARRAY_DOUBLE2(Density_Matrix[block], Info->Dim_Block[block]);
      FREE_ARRAY_DOUBLE2(Vector_DM[block]     , Info->Dim_Block[block]);
      FREE_ARRAY_DOUBLE1(Value_DM[block]);
   }
   
   free(Density_Matrix);
   free(Vector_DM);
   free(Value_DM);
   
   FREE_ARRAY_INT1(Block_Sorted);
   FREE_ARRAY_INT1(Basis_Sorted);
   FREE_ARRAY_INT1(Q_Number_Sorted);
   FREE_ARRAY_INT1(Row_Base);
   
   DMRG_FREE_DENSITY_MATRIX_INFO(Info);
   
   Dmrg_Time->make_dens_mat = omp_get_wtime() - Dmrg_Time->make_dens_mat;
   
   return Dmrg_System;
   
}
