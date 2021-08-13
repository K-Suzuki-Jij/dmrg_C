#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "dmrg.h"

void DMRG_EXTEND_AND_REDUCE(CRS1 *M_LL, CRS1 *M_On, CRS1 *Out, int *Ele_LL, char Type[], char Sign[], CCS1 *T_M, CRS1 *T_MD, CRS1 *Work_CRS, CCS1 *Work_CCS, double *Work_Vec, DMRG_BASIS_LLLR *Dmrg_Basis) {
   
   if(strcmp(Type, "LL") == 0){
      DMRG_MAKE_LL_OP_LLLR(M_LL, Dmrg_Basis, Work_CRS);
      CRS_CRS_CCS_PRODUCT(T_MD, Work_CRS, T_M, Out, Work_CCS, Work_Vec);
   }
   
   else if (strcmp(Type, "LR") == 0) {
      DMRG_MAKE_LR_OP_LLLR(M_On, Ele_LL, Sign, Dmrg_Basis, Work_CRS);
      CRS_CRS_CCS_PRODUCT(T_MD, Work_CRS, T_M, Out, Work_CCS, Work_Vec);
   }
   
   else if (strcmp(Type, "LL_LR") == 0) {
      DMRG_MAKE_LLLR_OP_LLLR(M_LL, M_On, 1.0, Ele_LL, Sign, Type, Dmrg_Basis, Work_CRS);
      CRS_CRS_CCS_PRODUCT(T_MD, Work_CRS, T_M, Out, Work_CCS, Work_Vec);
   }
   
   else {
      printf("Error in DMRG_EXTEND_AND_REDUCE\n");
      printf("Type = %s\n", Type);
      exit(1);
   }
   
}
