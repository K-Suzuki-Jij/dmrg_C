#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <lapack.h>
#include <sys/stat.h>
#include "SML.h"

void LANCZOS_SYM(BOX_LAN *Box) {
   
   //Check input matrix
   if (Box->M->row_dim != Box->M->col_dim) {
      printf("Error in LANCZOS_SYM\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->M->row_dim <= 0 || Box->M->col_dim <= 0) {
      printf("Error in LANCZOS_SYM\n");
      printf("row_dim(%d) or col_dim(%d) is illegal value\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->max_step <= Box->min_step) {
      printf("Error in LANCZOS_SYM\n");
      printf("Box->max_step=%d,Box->min_step=%d\n", Box->max_step, Box->min_step);
      exit(1);
   }
   
   long i,j;
   int m_size = Box->M->row_dim;
   double temp;
   FILE *file;
   
   //Diagonalize by lapack if m_size < 1000
   if (m_size < 1000) {
      double **Temp_Eigen_Vector = GET_ARRAY_DOUBLE2(1, m_size);
      LAPACK_DSYEV_CRS1(Box->M, Box->eig_val, Temp_Eigen_Vector, 1, 1);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/LANCZOS_SYM_Step.txt","a+")) == NULL) {
         printf("Error in LANCZOS_SYM\n");
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "%d   Caluclated by LAPACK\n", 0);
      fclose(file);
      
      for(i = 0; i < m_size; i++) {
         Box->Eig_Vec[i] = Temp_Eigen_Vector[0][i];
      }
      FREE_ARRAY_DOUBLE2(Temp_Eigen_Vector, 1);
      
      return;
   }
   
   //LANCZOS_SYM
   double *Vector_0     = GET_ARRAY_DOUBLE1(m_size);
   double *Vector_1     = GET_ARRAY_DOUBLE1(m_size);
   double *Vector_2     = GET_ARRAY_DOUBLE1(m_size);
   double **Temp_V      = GET_ARRAY_DOUBLE2(Box->p_threads, m_size);
   double **Rits_Vector = GET_ARRAY_DOUBLE2(Box->max_block_step+1,m_size);
   double *Diag         = GET_ARRAY_DOUBLE1(Box->max_block_step);
   double *Off_Diag     = GET_ARRAY_DOUBLE1(Box->max_block_step);
   
   //LAPACK
   double *Lap_D = GET_ARRAY_DOUBLE1(Box->max_block_step+1);
   double *Lap_E = GET_ARRAY_DOUBLE1(Box->max_block_step+1);
   double *Lap_V = GET_ARRAY_DOUBLE1((Box->max_block_step+1)*(Box->max_block_step+1));
   double *Work  = GET_ARRAY_DOUBLE1((Box->max_block_step+1)*2);
   int n = 0,ldz = 0,info = 0;
   char jobz = 'V';
   
   //Convergence check
   int conv_flag;
   double E = pow(10,15);
   
   if (strcmp(Box->Guess, "No") == 0) {
      //Set initial a vector
      srand((unsigned int)time(NULL));
#pragma omp parallel for num_threads (Box->p_threads)
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
   else {
      printf("Error in LANCZOS_SYM\n");
      printf("Box->Guess=%s is error\n", Box->Guess);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "Error in LANCZOS_SYM\n");
      fprintf(file, "Box->Guess=%s is error\n", Box->Guess);
      fclose(file);
      exit(1);
   }
   
   int total_step = 0;
   while(total_step < Box->max_step){
      
      //Normalize initial vectors
      NORMALIZE(Vector_0, m_size, Box->p_threads);
      
      COPY_DOUBLE1(Vector_0, Rits_Vector[0], m_size, Box->p_threads);
      
      MATRIX_VECTOR_PRODUCT_SYM(Box->M, Vector_0, Vector_1, Temp_V, Box->p_threads);
      
      Diag[0] = INNER_PRODUCT(Vector_0, Vector_1, m_size, Box->p_threads);
      
      temp = Diag[0];
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         Vector_1[i] = Vector_1[i] - temp*Vector_0[i];
      }
      
      //start iteration
      int block_step;
      for (block_step = 1; (block_step < Box->max_block_step) && (total_step < Box->max_step); block_step++) {
        
         total_step = total_step + 1;
         
         COPY_DOUBLE1(Vector_1, Vector_2, m_size, Box->p_threads);
         
         Off_Diag[block_step] = sqrt(INNER_PRODUCT(Vector_2, Vector_2, m_size, Box->p_threads));
         
         NORMALIZE(Vector_2, m_size, Box->p_threads);
         
         COPY_DOUBLE1(Vector_2, Rits_Vector[block_step], m_size, Box->p_threads);
         
         MATRIX_VECTOR_PRODUCT_SYM(Box->M, Vector_2, Vector_1, Temp_V, Box->p_threads);
         
         Diag[block_step] = INNER_PRODUCT(Vector_1, Vector_2, m_size, Box->p_threads);
         
         //Diagonalize by LAPACK
         if (total_step >= Box->min_step || block_step == Box->max_block_step - 1) {
            
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
               printf("\rLANCZOS_error[%d]=%e", total_step, fabs(Lap_V[n-1]));
               conv_flag = 0;
               if (fabs(Lap_V[n-1]) < Box->acc) {
                  conv_flag = 1;
               }
            }
            else if (strcmp(Box->Type, "Normal") == 0) {
               fflush(stdout);
               printf("\rLANCZOS_error[%d]=%e", total_step, fabs(E - Lap_D[0]));
               conv_flag = 0;
               if (fabs(E - Lap_D[0]) < Box->acc) {
                  conv_flag = 1;
               }
            }
            else {
               printf("Error in LANCZOS_SYM\n");
               printf("Box->Type=%s is error\n", Box->Type);
               mkdir("SML_out",0777);
               if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                  printf("Can't open file.\n");
                  exit(1);
               }
               fprintf(file, "Error in LANCZOS_SYM\n");
               fprintf(file, "Box->Type=%s is error\n", Box->Type);
               fclose(file);
               exit(1);
            }
            
            if (conv_flag == 1) {
               
               Box->eig_val[0] = Lap_D[0];
               
#pragma omp parallel for private (j,temp) num_threads (Box->p_threads)
               for (i = 0; i < m_size; i++) {
                  temp = 0;
                  for (j = 0; j < n; j++) {
                     temp = temp + Lap_V[j]*Rits_Vector[j][i];
                  }
                  Box->Eig_Vec[i] = temp;
               }
             
               //Normalize eigen vectors (just in case)
               NORMALIZE(Box->Eig_Vec, m_size, Box->p_threads);
               
               mkdir("SML_out",0777);
               if ((file = fopen("./SML_out/LANCZOS_SYM_Step.txt","a+")) == NULL) {
                  printf("Error in LANCZOS_SYM\n");
                  printf("Can't open file.\n");
                  exit(1);
               }
               fprintf(file, "%d\n", total_step);
               fclose(file);
               
               FREE_ARRAY_DOUBLE1(Vector_0);
               FREE_ARRAY_DOUBLE1(Vector_1);
               FREE_ARRAY_DOUBLE1(Vector_2);
               FREE_ARRAY_DOUBLE2(Rits_Vector, Box->max_block_step+1);
               FREE_ARRAY_DOUBLE1(Diag);
               FREE_ARRAY_DOUBLE1(Off_Diag);
               FREE_ARRAY_DOUBLE1(Lap_D);
               FREE_ARRAY_DOUBLE1(Lap_E);
               FREE_ARRAY_DOUBLE1(Lap_V);
               FREE_ARRAY_DOUBLE1(Work);
               FREE_ARRAY_DOUBLE2(Temp_V, Box->p_threads);
               
               printf("\r                           ");
               printf("                           \r");
               
               return;
            }
            
            if (block_step == Box->max_block_step - 1) {
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
      
#pragma omp parallel for private (j,temp) num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         temp = 0;
         for (j = 0; j < n; j++) {
            temp = temp + Lap_V[j]*Rits_Vector[j][i];
         }
         Vector_0[i] = temp;
      }
      
   }
   
   printf("Error in LANCZOS_SYM\n");
   printf("Does not converge (max step is %d)\n", Box->max_step);
   mkdir("SML_out",0777);
   if((file = fopen("./SML_out/error.txt","a+")) == NULL){
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "Error in LANCZOS_SYM\n");
   fprintf(file, "Does not converge (max step is %d)\n", Box->max_step);
   fclose(file);
   exit(1);
   
}







