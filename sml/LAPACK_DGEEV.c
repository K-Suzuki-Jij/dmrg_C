#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <lapack.h>
#include <math.h>
#include "SML.h"

void LAPACK_DGEEV(double **Ham,
                  int dim,
                  double *Eigen_Value_Real,
                  double *Eigen_Value_Img,
                  double **Eigen_Vector_Real,
                  double **Eigen_Vector_Img,
                  int eigval_num,
                  int eigvec_num
                  ){
   
   //If eigenvalues of Ham are complex, eigen vectors are not caluculated
   
   int n,lda,ldvl,ldvr,lwork,info;
   char jobvl,jobvr;
   double zero;
   
   jobvl = 'N';
   jobvr = 'V';
   n = dim;
   zero = pow(10,-16);
   
   double *A = GET_ARRAY_DOUBLE1(n*n);
   
   lda = n;
   
   double *WR = GET_ARRAY_DOUBLE1(n);
   double *WI = GET_ARRAY_DOUBLE1(n);
   double *VL = GET_ARRAY_DOUBLE1(n*n);
   double *VR = GET_ARRAY_DOUBLE1(n*n);
   
   ldvl = n;
   ldvr = n;
   
   double *WORK = GET_ARRAY_DOUBLE1(4*n);
   
   lwork = 4*n;
   
   int i,j;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         A[j*n + i] = Ham[i][j];
      }
   }
   
   dgeev_(&jobvl, &jobvr, &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, &info);
   
   
   for (i = 0; i < eigval_num; i++) {
      Eigen_Value_Real[i] = WR[i];
      Eigen_Value_Img[i] = WI[i];
   }
   
   i = 0;
   while (i < eigvec_num) {
      if (fabs(WI[i]) < zero) {
         for (j = 0; j < n; j++) {
            Eigen_Vector_Real[i][j] = VR[i*n + j];
            Eigen_Vector_Img[i][j] = 0;
         }
         i = i + 1;
      }
      else {
         for (j = 0; j < n; j++) {
            Eigen_Vector_Real[i][j] = VR[i*n + j];
            Eigen_Vector_Img[i][j] = VR[(i + 1)*n + j];
         }
         if (i + 1 == eigvec_num) {
            break;
         }
         for (j = 0; j < n; j++) {
            Eigen_Vector_Real[i + 1][j] = VR[i*n + j];
            Eigen_Vector_Img[i + 1][j] = -VR[(i + 1)*n + j];
         }
         i = i + 2;
      }
   }
   
   free(A);
   free(WR);
   free(WI);
   free(VL);
   free(VR);
   free(WORK);

}
