#include <stdio.h>
#include <stdlib.h>
//#include <lapack.h>
#include "SML.h"

extern void dgesvd_(char *, char *, int *, int *, double *, int *, double *, double *, int *, double *, int *, double *, int *, int *);

void LAPACK_DGESVD(double **M,
                   int row,
                   int col,
                   double *Eigen_Value,
                   double **Eigen_Vector_R,
                   double **Eigen_Vector_L,
                   int eigval_num,
                   int eigvec_num_r,
                   int eigvec_num_l
                   ){
   
   int i,j;
   int m,n,lda,ldu,ldvt,lwork,info;
   char jobu,jobvt;
   
   jobu = 'A';
   jobvt = 'A';
   m = row;
   n = col;
   
   double *A = GET_ARRAY_DOUBLE1(m*n);

   for (j = 0; j < n; j++) {
      for (i = 0; i < m; i++) {
         A[j*n + i] = M[i][j];
      }
   }
   
   lda = m;
   
   double *S = GET_ARRAY_DOUBLE1(m);
  
   double *U = GET_ARRAY_DOUBLE1(m*m);
   
   ldu = m;
   
   double *VT = GET_ARRAY_DOUBLE1(n*n);
   
   ldvt = n;
   
   lwork = 9*m;
   
   double *WORK = GET_ARRAY_DOUBLE1(lwork);
   
   dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, WORK, &lwork, &info);
   
   for (i = 0; i < eigval_num; i++) {
      Eigen_Value[i] = S[i];
   }
   
   for (i = 0; i < eigvec_num_r; i++) {
      for (j = 0; j < n; j++) {
         Eigen_Vector_R[i][j] = VT[i*n + j];
      }
   }
   
   
   for (i = 0; i < eigvec_num_l; i++) {
      for (j = 0; j < m; j++) {
         Eigen_Vector_L[i][j] = U[j*m + i];
      }
   }
   
   free(A);
   free(S);
   free(U);
   free(VT);
   free(WORK);
   
}
