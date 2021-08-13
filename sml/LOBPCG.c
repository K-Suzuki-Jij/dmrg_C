#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
//#include <lapack.h>
#include <time.h>
#include <sys/stat.h>
#include "SML.h"

extern void dspgv_(int *, char *, char *, int *, double *, double *, double *, double *, int *, double *, int *);

void LOBPCG(BOX_LOBPCG *Box) {
   
   //Check input matrix
   if (Box->M->row_dim != Box->M->col_dim) {
      printf("Error in LOBPCG\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->M->row_dim <= 0 || Box->M->col_dim <= 0) {
      printf("Error in LOBPCG\n");
      printf("row_dim(%d) or col_dim(%d) is illegal value\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   
   long i;
   int m_size = Box->M->row_dim;
   double alpha,beta;
   FILE *file;
   
   //Diagonalize by lapack if m_size < 1000
   if (Box->M->row_dim < 1000) {
      double **Temp_Eigen_Vector = GET_ARRAY_DOUBLE2(1, m_size);
      LAPACK_DSYEV_CRS1(Box->M, Box->eig_val, Temp_Eigen_Vector, 1, 1);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/LOBPCG-Step.txt","a+")) == NULL) {
         printf("Error in LOBPCG\n");
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "%d   Caluclated by LAPACK\n", 0);
      fclose(file);
      for (i = 0; i < m_size; i++) {
         Box->Eig_Vec[i] = Temp_Eigen_Vector[0][i];
      }
      FREE_ARRAY_DOUBLE2(Temp_Eigen_Vector, 1);
      return;
   }
   
   //LOBPCG
   double *V0 = GET_ARRAY_DOUBLE1(m_size);
   double *V1 = GET_ARRAY_DOUBLE1(m_size);
   double *V2 = GET_ARRAY_DOUBLE1(m_size);
   double *W0 = GET_ARRAY_DOUBLE1(m_size);
   double *W1 = GET_ARRAY_DOUBLE1(m_size);
   double *W2 = GET_ARRAY_DOUBLE1(m_size);
   double *R  = GET_ARRAY_DOUBLE1(m_size);
   double zero = pow(10,-15);
   double temp0,temp1,temp2,E;
   
   //LAPACK
   double *L_AP   = GET_ARRAY_DOUBLE1(6);
   double *L_BP   = GET_ARRAY_DOUBLE1(6);
   double *L_W    = GET_ARRAY_DOUBLE1(3);
   double *L_Z    = GET_ARRAY_DOUBLE1(9);
   double *L_Work = GET_ARRAY_DOUBLE1(9);
   int L_itype,L_n,L_ldz,L_info;
   char L_jobz,L_uplo;
   L_itype = 1;
   L_jobz = 'V';
   L_uplo = 'U';
   
   //Set initial vectors
   srand((unsigned int)time(NULL));
   
#pragma omp parallel for num_threads (Box->p_threads)
   for (i = 0; i < m_size; i++) {
      V0[i] = rand()%100000 - 50000;
   }
   NORMALIZE(V0, m_size, Box->p_threads);
   
   MATRIX_VECTOR_PRODUCT(Box->M, V0, W0, Box->p_threads);
   
   int block_step;
   
   for (block_step = 0; block_step < Box->max_step; block_step++) {
      
      L_AP[0] = INNER_PRODUCT(V0, W0, m_size, Box->p_threads);
      L_AP[1] = INNER_PRODUCT(V0, W1, m_size, Box->p_threads);
      L_AP[2] = INNER_PRODUCT(V1, W1, m_size, Box->p_threads);
      L_AP[3] = INNER_PRODUCT(V0, W2, m_size, Box->p_threads);
      L_AP[4] = INNER_PRODUCT(V1, W2, m_size, Box->p_threads);
      L_AP[5] = INNER_PRODUCT(V2, W2, m_size, Box->p_threads);
      
      L_BP[0] = 1.0;//INNER_PRODUCT(V0,V0,m_size,Box->p_threads)
      L_BP[1] = INNER_PRODUCT(V0, V1, m_size, Box->p_threads);
      L_BP[2] = 1.0;//INNER_PRODUCT(V1,V1,m_size,Box->p_threads);
      L_BP[3] = INNER_PRODUCT(V0, V2, m_size, Box->p_threads);
      L_BP[4] = INNER_PRODUCT(V1, V2, m_size, Box->p_threads);
      L_BP[5] = 1.0;//INNER_PRODUCT(V2,V2,m_size,Box->p_threads);

      if (block_step >= 2) {
         L_n = 3;
      }
      else if (block_step == 1) {
         L_n = 2;
      }
      else {
         L_n = 1;
      }
      L_ldz = L_n;

      dspgv_(&L_itype, &L_jobz, &L_uplo, &L_n, L_AP, L_BP, L_W, L_Z, &L_ldz, L_Work, &L_info);
      temp0 = L_Z[0];
      temp1 = L_Z[1];
      temp2 = L_Z[2];
      E = L_W[0];
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         R[i] = temp0*W0[i] + temp1*W1[i] + temp2*W2[i] - E*(temp0*V0[i] + temp1*V1[i] + temp2*V2[i]);
      }
      
      //Check convergence
      if (INNER_PRODUCT(R, R, m_size, Box->p_threads) < Box->acc) {
         Box->eig_val[0] = E;
         for (i = 0; i < m_size; i++) {
            Box->Eig_Vec[i] = temp0*W0[i] + temp1*W1[i] + temp2*W2[i];
         }
         NORMALIZE(Box->Eig_Vec, m_size, Box->p_threads);
         
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/LOBPCG-Step.txt","a+")) == NULL) {
            printf("Error in LOBPCG\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", block_step);
         fclose(file);
         
         FREE_ARRAY_DOUBLE1(V0);
         FREE_ARRAY_DOUBLE1(V1);
         FREE_ARRAY_DOUBLE1(V2);
         FREE_ARRAY_DOUBLE1(W0);
         FREE_ARRAY_DOUBLE1(W1);
         FREE_ARRAY_DOUBLE1(W2);
         FREE_ARRAY_DOUBLE1(R);
         FREE_ARRAY_DOUBLE1(L_AP);
         FREE_ARRAY_DOUBLE1(L_BP);
         FREE_ARRAY_DOUBLE1(L_W);
         FREE_ARRAY_DOUBLE1(L_Z);
         FREE_ARRAY_DOUBLE1(L_Work);
         
         return;
      }
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         V0[i] = temp0*V0[i] + temp1*V1[i] + temp2*V2[i];
      }
      alpha = 1.0/sqrt(INNER_PRODUCT(V0, V0, m_size, Box->p_threads));
      
      NORMALIZE(V0, m_size, Box->p_threads);
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         V2[i] = temp1*V1[i] + temp2*V2[i];
      }
      beta = INNER_PRODUCT(V2, V2, m_size, Box->p_threads);
      
      if (beta > zero) {
         NORMALIZE(V2, m_size, Box->p_threads);
         beta = 1.0/sqrt(beta);
      }
      else {
         beta = 0;
      }
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         V1[i] = R[i];
      }
      
      NORMALIZE(V1, m_size, Box->p_threads);
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         W0[i] = alpha*(temp0*W0[i] + temp1*W1[i] + temp2*W2[i]);
      }
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         W2[i] = beta*(temp1*W1[i] + temp2*W2[i]);
      }
      
      MATRIX_VECTOR_PRODUCT(Box->M, V1, W1, Box->p_threads);
   }
   
   printf("Error in LOBPCG\n");
   printf("Does not converge (max step is %d)\n", Box->max_step);
   mkdir("SML_out",0777);
   if((file = fopen("./SML_out/error.txt","a+")) == NULL){
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file,"Error in LOBPCG\n");
   fprintf(file,"Does not converge (max step is %d)\n", Box->max_step);
   fclose(file);
   exit(1);
   
}
