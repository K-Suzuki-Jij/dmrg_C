#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/stat.h>
#include "SML.h"

void COPY_CCS1(CCS1 *M, CCS1 *Copyed, int p_threads) {
   
   long i,j;
   FILE *file;
   
   if (Copyed->max_col < M->col_dim) {
      printf("error in COPY_CCS1\n");
      printf("Nedd more Copyed->col_dim at least %d\n", M->col_dim);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "error in COPY_CCS1\n");
      fprintf(file, "Nedd more Copyed->col_dim at least %d\n", M->col_dim);
      fclose(file);
      exit(1);
   }
   if (Copyed->max_val < M->Col[M->col_dim]) {
      printf("error in COPY_CCS1\n");
      printf("Nedd more Copyed->max_val at least %ld\n", M->Col[M->col_dim]);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file,"error in COPY_CCS1\n");
      fprintf(file,"Nedd more Copyed->max_val at least %ld\n", M->Col[M->col_dim]);
      fclose(file);
      exit(1);
   }
   
#pragma omp parallel for private (j) num_threads (p_threads)
   for (i = 0; i < M->col_dim; i++) {
      for (j = M->Col[i]; j < M->Col[i+1]; j++) {
         Copyed->Row[j] = M->Row[j];
         Copyed->Val[j] = M->Val[j];
      }
      Copyed->Col[i+1] = M->Col[i+1];
   }
   
   Copyed->row_dim = M->row_dim;
   Copyed->col_dim = M->col_dim;
   
}
