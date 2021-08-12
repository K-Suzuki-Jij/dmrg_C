#include <stdio.h>
#include <stdlib.h>
#include <lapack.h>
#include "SML.h"

void LAPACK_DSYEV_CRS1(CRS1 *Ham,
                       double *Eigen_Value,
                       double **Eigen_Vector,
                       int eigval_num,
                       int eigvec_num
                       ){
   
   long i,j;
   int m_size = Ham->row_dim;
   char jobz='V';
   char uplo='U';
   int n,lda,lwork,info;
   
   double *Mat   = GET_ARRAY_DOUBLE1(m_size*m_size);
   double *Value = GET_ARRAY_DOUBLE1(m_size);
   double *Work  = GET_ARRAY_DOUBLE1(3*m_size);
   
   for (i = 0; i < m_size; i++) {
      for (j = Ham->Row[i]; j < Ham->Row[i+1]; j++) {
         Mat[i*m_size + Ham->Col[j]] = Ham->Val[j];
         Mat[Ham->Col[j]*m_size + i] = Ham->Val[j]; //
      }
   }
   
   /*
   for (i = 0; i < m_size; i++) {
      for (j = 0; j < m_size; j++) {
         if (fabs(Mat[i*m_size + j]) > pow(10,-15)) {
         printf("%-2ld %-2ld  %-.15lf\n",i,j,Mat[i*m_size + j]);
         }
         
      }
   }
   
   exit(1);
    */
   
   n = m_size;
   lda = m_size;
   lwork = 3*m_size;
   
   dsyev_(&jobz, &uplo, &n, Mat, &lda, Value, Work, &lwork, &info);
   
   for (i = 0; i < eigvec_num; i++) {
      for (j = 0; j < m_size; j++) {
         Eigen_Vector[i][j] = Mat[j+m_size*i];
      }
   }
   
   for (i = 0; i < eigval_num; i++) {
      Eigen_Value[i] = Value[i];
   }
   
   FREE_ARRAY_DOUBLE1(Mat);
   FREE_ARRAY_DOUBLE1(Value);
   FREE_ARRAY_DOUBLE1(Work);
   
}
