#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "dmrg.h"

double **DMRG_GET_DENSITY_MATRIX(int block, int *Dim_Block, DMRG_BASIS *Dmrg_Basis, double *Vec, int p_threads) {
   
   int dim_matrix = Dim_Block[block];
   int dim_LLLR   = Dmrg_Basis->dim_LLLR;
   int dim_RR     = Dmrg_Basis->dim_RR;
   int dim_onsite = Dmrg_Basis->dim_onsite;
   
   double **M = GET_ARRAY_DOUBLE2(dim_matrix, dim_matrix);
   int row,col,row_base,col_base,iter,base,sum,inv_row,inv_col;
   int LL_row,LR_row,LL_col,LR_col,RR_row,RL_row,RR_col,RL_col;
   int *Base_LLLR = GET_ARRAY_INT1(dim_LLLR + 1);
   long inv_sup_row,inv_sup_col;
   
   for (iter = 0; iter < dim_LLLR; iter++) {
      Base_LLLR[iter + 1] = Base_LLLR[iter] + Dmrg_Basis->Sum_LLLR[iter];
   }
   
   //Check Point
   if (Base_LLLR[dim_LLLR] != Dmrg_Basis->dim_LLLRRRRL) {
      printf("Error in DMRG_GET_DENSITY_MATRIX\n");
      printf("Base_LLLR[%d] (%d) != dim_LLLRRRRL (%d)\n", dim_LLLR, Base_LLLR[dim_LLLR], Dmrg_Basis->dim_LLLRRRRL);
      exit(1);
   }
   
   base = 0;
   for (iter = 0; iter < block; iter++) {
      base = base + Dim_Block[iter];
   }
   
   
#pragma omp parallel for private (row_base,LL_row,LR_row,col,col_base,LL_col,LR_col,sum,iter,RR_row,RR_col,RL_row,RL_col,inv_row,inv_col,inv_sup_row,inv_sup_col) num_threads (p_threads)
   for (row = 0; row < dim_matrix; row++) {
      row_base = row + base;
      LL_row = Dmrg_Basis->LL_LLLR[row_base];
      LR_row = Dmrg_Basis->LR_LLLR[row_base];
      for (col = row; col < dim_matrix; col++) {
         col_base = col + base;
         
         //Check Point
         if (row_base >= dim_LLLR || col_base >= dim_LLLR) {
            printf("Error in DMRG_GET_DENSITY_MATRIX\n");
            printf("row_base (%d) >= dim_LLLR (%d) || col_base (%d) >= dim_LLLR (%d)\n",
                   row_base, dim_LLLR, col_base, dim_LLLR);
            exit(1);
         }
         
         //Check Point
         if (Dmrg_Basis->Sum_LLLR[row_base] != Dmrg_Basis->Sum_LLLR[col_base]) {
            printf("Error in DMRG_GET_DENSITY_MATRIX\n");
            printf("Sum_LLLR[%d] (%d) != Sum_LLLR[%d] (%d)\n",
                   row_base, Dmrg_Basis->Sum_LLLR[row_base], col_base, Dmrg_Basis->Sum_LLLR[col_base]);
            exit(1);
         }
         
         LL_col = Dmrg_Basis->LL_LLLR[col_base];
         LR_col = Dmrg_Basis->LR_LLLR[col_base];
         sum = Dmrg_Basis->Sum_LLLR[row_base];
         
         for (iter = 0; iter < sum; iter++) {
            RR_row = Dmrg_Basis->RR_LLLRRRRL[Base_LLLR[row_base] + iter];
            RR_col = Dmrg_Basis->RR_LLLRRRRL[Base_LLLR[col_base] + iter];
            RL_row = Dmrg_Basis->RL_LLLRRRRL[Base_LLLR[row_base] + iter];
            RL_col = Dmrg_Basis->RL_LLLRRRRL[Base_LLLR[col_base] + iter];
            
            //Check Point
            if (RR_row != RR_col || RL_row != RL_col) {
               printf("Error in DMRG_GET_DENSITY_MATRIX\n");
               printf("RR_row (%d) != RR_col (%d) || RL_row (%d) != RL_col (%d)\n",
                      RR_row, RR_col, RL_row, RL_col
                      );
               exit(1);
            }
            inv_sup_row = (long)LL_row*dim_onsite*dim_RR*dim_onsite + (long)LR_row*dim_RR*dim_onsite + RR_row*dim_onsite + RL_row;
            inv_sup_col = (long)LL_col*dim_onsite*dim_RR*dim_onsite + (long)LR_col*dim_RR*dim_onsite + RR_col*dim_onsite + RL_col;
            inv_row = Dmrg_Basis->Inv_LLLRRRRL[inv_sup_row];
            inv_col = Dmrg_Basis->Inv_LLLRRRRL[inv_sup_col];
            
            //Check Point
            if (inv_row < 0 || inv_col < 0) {
               printf("Error in DMRG_GET_DENSITY_MATRIX\n");
               printf("inv_row (%d) < 0 || inv_col (%d) < 0\n",
                      inv_row, inv_col
                      );
               exit(1);
            }
            
            M[row][col] = M[row][col] + Vec[inv_row]*Vec[inv_col];
         }
      }
   }
   
   for (row = 0; row < dim_matrix; row++) {
      for (col = row; col < dim_matrix; col++) {
         M[col][row] = M[row][col];
      }
   }
   
   FREE_ARRAY_INT1(Base_LLLR);

   return M;
   
}
