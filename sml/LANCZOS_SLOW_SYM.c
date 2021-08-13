#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <sys/stat.h>
//#include <lapack.h>
#include "SML.h"

extern void dstev_(char *, int *, double *, double *, double *, int *, double *, int *);

void LANCZOS_SLOW_SYM(BOX_LAN *Box) {
   
   //Check input matrix
   if (Box->M->row_dim != Box->M->col_dim) {
      printf("Error in LANCZOS_SLOW_SYM\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->M->row_dim <= 0 || Box->M->col_dim <= 0) {
      printf("Error in LANCZOS_SLOW_SYM\n");
      printf("row_dim(%d) or col_dim(%d) is illegal value\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->max_step <= Box->min_step) {
      printf("Error in LANCZOS_SLOW_SYM\n");
      printf("Box->max_step=%d,Box->min_step=%d\n", Box->max_step, Box->min_step);
      exit(1);
   }
   
   long i;
   int m_size = Box->M->row_dim;
   double temp;
   FILE *file;
   
   //Diagonalize by lapack if m_size < 1000
   if (m_size < 1000) {
      double **Temp_Eigen_Vector = GET_ARRAY_DOUBLE2(1, m_size);
      LAPACK_DSYEV_CRS1(Box->M, Box->eig_val, Temp_Eigen_Vector, 1, 1);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/LANCZOS_SLOW_SYM_Step.txt","a+")) == NULL) {
         printf("Error in LANCZOS_SLOW_SYM\n");
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
   
   //LANCZOS_SLOW_SYM
   double *Vector_0 = GET_ARRAY_DOUBLE1(m_size);
   double *Vector_1 = GET_ARRAY_DOUBLE1(m_size);
   double *Vector_2 = GET_ARRAY_DOUBLE1(m_size);
   double **Temp_V  = GET_ARRAY_DOUBLE2(Box->p_threads, m_size);
   double *Diag     = GET_ARRAY_DOUBLE1(Box->max_step);
   double *Off_Diag = GET_ARRAY_DOUBLE1(Box->max_step);
   unsigned int seed = (unsigned int)time(NULL);
   
   //LAPACK
   double *Lap_D = GET_ARRAY_DOUBLE1(Box->max_step+1);
   double *Lap_E = GET_ARRAY_DOUBLE1(Box->max_step+1);
   double *Lap_V = GET_ARRAY_DOUBLE1((Box->max_step+1)*(Box->max_step+1));
   double *Work  = GET_ARRAY_DOUBLE1((Box->max_step+1)*2);
   int n,ldz,info;
   char jobz = 'V';
   
   //Convergence check
   int conv_flag = 0;
   double E = pow(10,15);
   
   if (strcmp(Box->Guess, "No") == 0) {
      //Set initial a vector
      srand(seed);
      for (i = 0; i < m_size; i++) {
         Vector_0[i] = rand()%100000000 - 50000000;
      }
   }
   else if(strcmp(Box->Guess, "Yes") == 0){
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         Vector_0[i] = Box->Eig_Vec[i];
      }
   }
   else {
      printf("Error in LANCZOS_SLOW_SYM\n");
      printf("Box->Guess=%s is error\n",Box->Guess);
      exit(1);
   }
   
   //Normalize initial vectors
   NORMALIZE(Vector_0, m_size, Box->p_threads);
   MATRIX_VECTOR_PRODUCT_SYM(Box->M, Vector_0, Vector_1, Temp_V, Box->p_threads);
   Diag[0] = INNER_PRODUCT(Vector_0, Vector_1, m_size, Box->p_threads);
   
   temp = Diag[0];
#pragma omp parallel for num_threads (Box->p_threads)
   for (i = 0; i < m_size; i++) {
      Vector_1[i] = Vector_1[i] - temp*Vector_0[i];
   }
   
   
   
   //start iteration
   int block_step;
   for (block_step = 1; block_step < Box->max_step; block_step++) {
      
      COPY_DOUBLE1(Vector_1, Vector_2, m_size, Box->p_threads);
      
      Off_Diag[block_step] = sqrt(INNER_PRODUCT(Vector_2, Vector_2, m_size, Box->p_threads));
      
      NORMALIZE(Vector_2, m_size, Box->p_threads);
      
      MATRIX_VECTOR_PRODUCT_SYM(Box->M, Vector_2, Vector_1, Temp_V, Box->p_threads);
      
      Diag[block_step] = INNER_PRODUCT(Vector_1, Vector_2, m_size, Box->p_threads);
      
      //Diagonalize by LAPACK
      if (block_step >= Box->min_step) {
         
         n = block_step + 1;
         
         for (i = 0; i < n; i++) {
            Lap_D[i] = Diag[i];
         }
         for (i = 1; i < n; i++) {
            Lap_E[i-1] = Off_Diag[i];
         }
         
         ldz = n;
         dstev_(&jobz, &n, Lap_D, Lap_E, Lap_V, &ldz, Work, &info);
         
         //Check convergence
         if (strcmp(Box->Type, "Tight") == 0) {
            fflush(stdout);
            printf("\rLANCZOS_SLOW_error[%d]=%e", block_step, fabs(Lap_V[n-1]));
            if (fabs(Lap_V[n-1]) < Box->acc) {
               conv_flag = 1;
            }
         }
         else if (strcmp(Box->Type, "Normal") == 0) {
            fflush(stdout);
            printf("\rLANCZOS_SLOW_error[%d]=%e", block_step, fabs(E - Lap_D[0]));
            if (fabs(E - Lap_D[0]) < Box->acc) {
               conv_flag = 1;
            }
         }
         else {
            printf("Error in LANCZOS_SLOW_SYM\n");
            printf("Box->Type=%s is error\n", Box->Type);
            exit(1);
         }
         
         if(conv_flag == 1){
            
            Box->eig_val[0] = Lap_D[0];
            
            mkdir("SML_out",0777);
            if ((file = fopen("./SML_out/LANCZOS_SLOW_SYM_Step.txt","a+")) == NULL) {
               printf("Error in LANCZOS_SLOW_SYM\n");
               printf("Can't open file.\n");
               exit(1);
            }
            fprintf(file,  "%d\n", block_step);
            fclose(file);
            
            printf("\r                           ");
            printf("                           \r");
            
            break;
         }
         
         E = Lap_D[0];
         
      }
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         Vector_1[i] = Vector_1[i] - Diag[block_step]*Vector_2[i] - Off_Diag[block_step]*Vector_0[i];
      }
      
      COPY_DOUBLE1(Vector_2, Vector_0, m_size, Box->p_threads);
      
   }
   
   //Repeat Lanczos step
   if (conv_flag == 1) {
      
      if (strcmp(Box->Guess, "No") == 0) {
         //Set initial a vector
         srand(seed);
         for (i = 0; i < m_size; i++) {
            Vector_0[i] = rand()%100000000 - 50000000;
            Box->Eig_Vec[i] = 0;
         }
      }
      else if (strcmp(Box->Guess, "Yes") == 0) {
#pragma omp parallel for num_threads (Box->p_threads)
         for (i = 0; i < m_size; i++) {
            Vector_0[i] = Box->Eig_Vec[i];
            Box->Eig_Vec[i] = 0;
         }
      }

      //Normalize initial vectors
      NORMALIZE(Vector_0, m_size, Box->p_threads);
      
      temp = Lap_V[0];
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         Box->Eig_Vec[i] = Box->Eig_Vec[i] + temp*Vector_0[i];
      }
      
      MATRIX_VECTOR_PRODUCT_SYM(Box->M, Vector_0, Vector_1, Temp_V, Box->p_threads);
      
      temp = Diag[0];
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         Vector_1[i] = Vector_1[i] - temp*Vector_0[i];
      }
      
      int iter;
      
      for (iter = 1; iter <= block_step; iter++) {
         
         COPY_DOUBLE1(Vector_1, Vector_2, m_size, Box->p_threads);
         
         NORMALIZE(Vector_2, m_size, Box->p_threads);
         
         temp = Lap_V[iter];
#pragma omp parallel for num_threads (Box->p_threads)
         for (i = 0; i < m_size; i++) {
            Box->Eig_Vec[i] = Box->Eig_Vec[i] + temp*Vector_2[i];
         }
         fflush(stdout);
         printf("\rLANCZOS_SLOW_CALC_VEC=%d/%d", iter, block_step);
         
         if (iter == block_step) {
            printf("\r                           ");
            printf("                           \r");
            break;
         }
         MATRIX_VECTOR_PRODUCT_SYM(Box->M, Vector_2, Vector_1, Temp_V, Box->p_threads);

#pragma omp parallel for num_threads (Box->p_threads)
         for (i = 0; i < m_size; i++) {
            Vector_1[i] = Vector_1[i] - Diag[iter]*Vector_2[i] - Off_Diag[iter]*Vector_0[i];
         }
         
         COPY_DOUBLE1(Vector_2, Vector_0, m_size, Box->p_threads);
         
      }
      
      //Normalize eigen vectors
      NORMALIZE(Box->Eig_Vec, m_size, Box->p_threads);
      FREE_ARRAY_DOUBLE1(Vector_0);
      FREE_ARRAY_DOUBLE1(Vector_1);
      FREE_ARRAY_DOUBLE1(Vector_2);
      FREE_ARRAY_DOUBLE1(Diag);
      FREE_ARRAY_DOUBLE1(Off_Diag);
      FREE_ARRAY_DOUBLE1(Lap_D);
      FREE_ARRAY_DOUBLE1(Lap_E);
      FREE_ARRAY_DOUBLE1(Lap_V);
      FREE_ARRAY_DOUBLE1(Work);
      FREE_ARRAY_DOUBLE2(Temp_V, Box->p_threads);
      return;
      
   }
   
   printf("Error in LANCZOS_SLOW_SYM\n");
   printf("Does not converge (max step is %d)\n", Box->max_step);
   mkdir("SML_out",0777);
   if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file,"Error in LANCZOS_SLOW_SYM\n");
   fprintf(file,"Does not converge (max step is %d)\n", Box->max_step);
   fclose(file);
   exit(1);
   
}
