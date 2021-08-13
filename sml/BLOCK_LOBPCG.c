#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
//#include <lapack.h>
#include "SML.h"

extern void dspgv_(int *, char *, char *, int *, double *, double *, double *, double *, int *, double *, int *);

void BLOCK_LOBPCG(BOX_BLOCK_LOBPCG *Box){
   
   //Check input matrix
   if (Box->M->row_dim != Box->M->col_dim) {
      printf("Error in BLOCK_LOBPCG\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->M->row_dim <= 0 || Box->M->col_dim <= 0) {
      printf("Error in BLOCK_LOBPCG\n");
      printf("row_dim(%d) or col_dim(%d) is illegal value\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->conv_eig_num < Box->eig_num) {
      printf("Error in BLOCK_LOBPCG\n");
      printf("Box->conv_eig_num(%d) < Box->eig_num(%d)\n", Box->conv_eig_num, Box->eig_num);
      exit(1);
   }
   
   long i,j,k;
   int int_temp,block_step;
   int m_size = Box->M->row_dim;
   FILE *file;
   
   //Diagonalize by lapack if m_size < 1000
   if (Box->M->row_dim < 1000) {
      LAPACK_DSYEV_CRS1(Box->M, Box->Eig_Val, Box->Eig_Vec, Box->eig_num, Box->eig_num);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/BLOCK_LOBPCG-Step.txt","a+")) == NULL) {
         printf("Error in BLOCK_LOBPCG\n");
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file,"%d   Caluclated by LAPACK\n",0);
      fclose(file);
      return;
   }
   
   //BLOCK_LOBPCG
   double **V0     = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double **V1     = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double **V2     = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double **W0     = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double **W1     = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double **W2     = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double **Temp_V = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double **R      = GET_ARRAY_DOUBLE2(Box->conv_eig_num,m_size);
   double *Alpha   = GET_ARRAY_DOUBLE1(Box->conv_eig_num);
   double *Beta    = GET_ARRAY_DOUBLE1(Box->conv_eig_num);
   double zero     = pow(10,-15);
   
   //LAPACK
   int L_itype,L_n,L_ldz,L_info;
   char L_jobz,L_uplo;
   L_itype = 1;
   L_jobz = 'V';
   L_uplo = 'U';
   L_n = Box->conv_eig_num*3;
   
   double *L_AP   = GET_ARRAY_DOUBLE1((L_n*(L_n+1))/2);
   double *L_BP   = GET_ARRAY_DOUBLE1((L_n*(L_n+1))/2);
   double *L_W    = GET_ARRAY_DOUBLE1(L_n);
   double *L_Z    = GET_ARRAY_DOUBLE1(L_n*L_n);
   double *L_Work = GET_ARRAY_DOUBLE1(L_n*L_n);
   
   //check convergence
   double *Norm = GET_ARRAY_DOUBLE1(Box->conv_eig_num);
   int conv_flag;
   int reorth_check;
   
   //Set initial vectors
   srand((unsigned int)time(NULL));
   
   for (i = 0;i < Box->conv_eig_num ;i++) {
#pragma omp parallel for num_threads (Box->p_threads)
      for(j = 0; j < m_size; j++){
         V0[i][j] = rand()%10000000 - 5000000;
      }
   }
   
   ORTHOGONALIZATION(V0, Box->conv_eig_num, m_size, 0, Box->p_threads);
   
   for (i = 0; i < Box->conv_eig_num; i++) {
      MATRIX_VECTOR_PRODUCT(Box->M, V0[i], W0[i], Box->p_threads);
   }
   

   for (block_step = 0 ;block_step < Box->max_step; block_step++){
      
      //Reorthogonalization every 10 steps
      reorth_check = 0;
      if (block_step%10 == 0) {
         reorth_check = 1;
         ORTHOGONALIZATION(V0, Box->conv_eig_num, m_size, 0, Box->p_threads);
         for (i = 0; i < Box->conv_eig_num; i++) {
#pragma omp parallel for num_threads (Box->p_threads)
            for (j = 0; j < m_size; j++) {
               V2[i][j] = 0;
               W2[i][j] = 0;
            }
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j <= i; j++) {
            L_AP[j+(i*(i+1))/2] = INNER_PRODUCT(V0[j], W0[i], m_size, Box->p_threads);
            if(i == j){
               L_BP[j+(i*(i+1))/2] = 1;
            }
            else {
               L_BP[j+(i*(i+1))/2] = 0;
            }
         }
      }
      
      int_temp = (Box->conv_eig_num*(Box->conv_eig_num + 1))/2;
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j < Box->conv_eig_num; j++) {
            L_AP[int_temp + i*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V0[j], W1[i], m_size, Box->p_threads);
            L_BP[int_temp + i*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V0[j], V1[i], m_size, Box->p_threads);
         }
      }
      
      int_temp = (Box->conv_eig_num*(Box->conv_eig_num+1))/2 + Box->conv_eig_num;
      
      for(i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j <= i; j++) {
            L_AP[int_temp + i*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V1[j], W1[i], m_size, Box->p_threads);
            if(i == j){
               L_BP[int_temp + i*Box->conv_eig_num + j + (i*(i+1))/2] = 1;
            }
            else{
               L_BP[int_temp + i*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V1[j], V1[i], m_size, Box->p_threads);
            }
         }
      }
      
      int_temp = (Box->conv_eig_num*(Box->conv_eig_num+1)) + Box->conv_eig_num*Box->conv_eig_num;
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j < Box->conv_eig_num; j++) {
            L_AP[int_temp + i*2*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V0[j], W2[i], m_size, Box->p_threads);
            L_BP[int_temp + i*2*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V0[j], V2[i], m_size, Box->p_threads);
         }
      }
      
      int_temp = Box->conv_eig_num*(Box->conv_eig_num+1) + Box->conv_eig_num*(Box->conv_eig_num+1);
      
      for (i = 0; i < Box->conv_eig_num; i++){
         for(j = 0; j < Box->conv_eig_num; j++){
            L_AP[int_temp + i*2*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V1[j], W2[i], m_size, Box->p_threads);
            L_BP[int_temp + i*2*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V1[j], V2[i], m_size, Box->p_threads);
         }
      }
      
      int_temp = Box->conv_eig_num*(Box->conv_eig_num+1) + Box->conv_eig_num*(Box->conv_eig_num+2);
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j <= i; j++) {
            L_AP[int_temp + i*2*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V2[j], W2[i], m_size, Box->p_threads);
            if(i == j){
               L_BP[int_temp + i*2*Box->conv_eig_num + j + (i*(i+1))/2] = 1;
            }
            else {
               L_BP[int_temp + i*2*Box->conv_eig_num + j + (i*(i+1))/2] = INNER_PRODUCT(V2[j], V2[i], m_size, Box->p_threads);
            }
         }
      }
      
      if (block_step >= 2) {
         L_n = 3*Box->conv_eig_num;
      }
      else if (block_step == 1) {
         L_n = 2*Box->conv_eig_num;
      }
      else {
         L_n = 1*Box->conv_eig_num;
      }
      
      if (block_step >= 2 && reorth_check == 1) {
         L_n = 2*Box->conv_eig_num;
      }
      
      L_ldz = L_n;
      
      dspgv_(&L_itype, &L_jobz, &L_uplo, &L_n, L_AP, L_BP, L_W, L_Z, &L_ldz, L_Work, &L_info);
      
      for (i = 0; i < Box->conv_eig_num; i++) {
#pragma omp parallel for num_threads (Box->p_threads)
         for (j = 0; j < m_size; j++){
            R[i][j] = 0;
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j < Box->conv_eig_num; j++) {
#pragma omp parallel for num_threads (Box->p_threads)
            for (k = 0; k < m_size; k++) {
               R[i][k] = R[i][k] + L_Z[j + i*L_n]*W0[j][k] + L_Z[j+Box->conv_eig_num + i*L_n]*W1[j][k] + L_Z[j+2*Box->conv_eig_num + i*L_n]*W2[j][k] -
               L_W[i]*(L_Z[j + i*L_n]*V0[j][k] + L_Z[j+Box->conv_eig_num + i*L_n]*V1[j][k] + L_Z[j+2*Box->conv_eig_num + i*L_n]*V2[j][k]);
            }
         }
      }
      
      conv_flag = 0;
      for (i = 0; i < Box->conv_eig_num; i++) {
         Norm[i] = INNER_PRODUCT(R[i], R[i], m_size, Box->p_threads);
         if (Norm[i] < Box->acc) {
            conv_flag = conv_flag + 1;
         }
      }
      
      fflush(stdout);
      printf("\rBLOCK_LOBPCG_error[%d]=%e", block_step, Norm[Box->eig_num-1]);
      
      //Check convergence
      if (conv_flag >= Box->eig_num) {
         for (i = 0; i < Box->eig_num; i++) {
#pragma omp parallel for num_threads (Box->p_threads)
            for (j = 0; j < m_size; j++) {
               Box->Eig_Vec[i][j] = 0;
            }
         }
         for (i = 0; i < Box->eig_num; i++) {
            for (j = 0; j < Box->eig_num; j++) {
#pragma omp parallel for num_threads (Box->p_threads)
               for (k = 0; k < m_size; k++) {
                  Box->Eig_Vec[i][k] = Box->Eig_Vec[i][k] + L_Z[j + i*L_n]*W0[j][k] +
                  L_Z[j+Box->conv_eig_num + i*L_n]*W1[j][k] +
                  L_Z[j+2*Box->conv_eig_num + i*L_n]*W2[j][k];
               }
            }
         }
         
         for (i = 0; i < Box->eig_num; i++) {
            NORMALIZE(Box->Eig_Vec[i], m_size, Box->p_threads);
            Box->Eig_Val[i] = L_W[i];
         }
         
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/BLOCK_LOBPCG-Step.txt","a+")) == NULL) {
            printf("Error in BLOCK_LOBPCG\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", block_step);
         fclose(file);
         
         FREE_ARRAY_DOUBLE2(V0,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE2(V1,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE2(V2,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE2(W0,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE2(W1,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE2(W2,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE2(Temp_V,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE2(R,Box->conv_eig_num);
         FREE_ARRAY_DOUBLE1(Alpha);
         FREE_ARRAY_DOUBLE1(Beta);
         
         FREE_ARRAY_DOUBLE1(L_AP);
         FREE_ARRAY_DOUBLE1(L_BP);
         FREE_ARRAY_DOUBLE1(L_W);
         FREE_ARRAY_DOUBLE1(L_Z);
         FREE_ARRAY_DOUBLE1(L_Work);
         
         FREE_ARRAY_DOUBLE1(Norm);
         printf("\r                        ");
         printf("                        \r");
         return;
      }
      
      for (i = 0;i < Box->conv_eig_num;i++) {
#pragma omp parallel for num_threads (Box->p_threads)
         for (j = 0; j < m_size; j++) {
            Temp_V[i][j] = V0[i][j];
            V0[i][j] = 0;
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j < Box->conv_eig_num; j++) {
#pragma omp parallel for num_threads (Box->p_threads)
            for (k = 0; k < m_size; k++) {
               V0[i][k] = V0[i][k] + L_Z[j + i*L_n]*Temp_V[j][k] + L_Z[j+Box->conv_eig_num + i*L_n]*V1[j][k] + L_Z[j+2*Box->conv_eig_num + i*L_n]*V2[j][k];
            }
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         Alpha[i] = 1.0/sqrt(INNER_PRODUCT(V0[i], V0[i], m_size, Box->p_threads));
         NORMALIZE(V0[i], m_size, Box->p_threads);
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
#pragma omp parallel for num_threads (Box->p_threads)
         for (j = 0; j < m_size; j++) {
            Temp_V[i][j] = V2[i][j];
            V2[i][j] = 0;
         }
      }
      
      for (i = 0;i < Box->conv_eig_num; i++) {
         for (j = 0;j < Box->conv_eig_num; j++) {
#pragma omp parallel for num_threads (Box->p_threads)
            for (k = 0; k < m_size; k++) {
               V2[i][k] = V2[i][k] +
               L_Z[j+Box->conv_eig_num + i*L_n]*V1[j][k] +
               L_Z[j+2*Box->conv_eig_num + i*L_n]*Temp_V[j][k];
            }
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         Beta[i] = INNER_PRODUCT(V2[i], V2[i], m_size, Box->p_threads);
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         if (Beta[i] > zero) {
            NORMALIZE(V2[i], m_size, Box->p_threads);
            Beta[i] = 1.0/sqrt(Beta[i]);
         }
         else {
            Beta[i] = 0;
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
#pragma omp parallel for num_threads (Box->p_threads)
         for (j = 0; j < m_size; j++) {
            V1[i][j] = R[i][j];
         }
         NORMALIZE(V1[i], m_size, Box->p_threads);
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
#pragma omp parallel for num_threads (Box->p_threads)
         for (j = 0; j < m_size; j++) {
            Temp_V[i][j] = W0[i][j];
            W0[i][j] = 0;
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j < Box->conv_eig_num; j++) {
#pragma omp parallel for num_threads (Box->p_threads)
            for (k = 0; k < m_size; k++) {
               W0[i][k] = W0[i][k] + Alpha[i]*(L_Z[j + i*L_n]*Temp_V[j][k] + L_Z[j+Box->conv_eig_num + i*L_n]*W1[j][k] + L_Z[j+2*Box->conv_eig_num + i*L_n]*W2[j][k]);
            }
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
#pragma omp parallel for num_threads (Box->p_threads)
         for (j = 0; j < m_size; j++) {
            Temp_V[i][j] = W2[i][j];
            W2[i][j] = 0;
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         for (j = 0; j < Box->conv_eig_num; j++) {
#pragma omp parallel for num_threads (Box->p_threads)
            for (k = 0; k < m_size; k++) {
               W2[i][k] = W2[i][k] + Beta[i]*(L_Z[j+Box->conv_eig_num + i*L_n]*W1[j][k] + L_Z[j+2*Box->conv_eig_num + i*L_n]*Temp_V[j][k]);
            }
         }
      }
      
      for (i = 0; i < Box->conv_eig_num; i++) {
         MATRIX_VECTOR_PRODUCT(Box->M, V1[i], W1[i], Box->p_threads);
      }
   }
   
   printf("Error in BLOCK_LOBPCG\n");
   printf("Does not converge (max step is %d)\n", Box->max_step);
   mkdir("SML_out",0777);
   if((file = fopen("./SML_out/error.txt","a+")) == NULL){
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file,"Error in BLOCK_LOBPCG\n");
   fprintf(file,"Does not converge (max step is %d)\n",Box->max_step);
   fclose(file);
   exit(1);
   
}
