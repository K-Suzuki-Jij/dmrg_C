#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "SML.h"
#include "dmrg.h"

double DMRG_CONVERGE_CHECK_OBC(double eig_val, CRS1 *M_LLLR, CRS1 *M_LRRL, CRS1 *M_LRRL_Sign, CRS1 *M_RRRL, double *Eig_V, double *Temp_V, int *Ele_RR, int p_threads, DMRG_BASIS *Dmrg_Basis) {
   
   DMRG_MATRIX_VECTOR_PRODUCT_OBC(M_LLLR, M_LRRL, M_LRRL_Sign, M_RRRL, Eig_V, Temp_V, Ele_RR, p_threads, Dmrg_Basis);
   
   int i;
   double norm = 0;
   
#pragma omp parallel for reduction (+:norm) num_threads (p_threads)
   for (i = 0; i < Dmrg_Basis->dim_LLLRRRRL; i++) {
      norm = norm + fabs(eig_val*Eig_V[i] - Temp_V[i]);
   }
   
   return norm;
}
