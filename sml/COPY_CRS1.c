#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/stat.h>
#include "SML.h"

void COPY_CRS1(CRS1 *M, CRS1 *Copyed, int p_threads) {
   
   long i,j;
   FILE *file;
   
   if (Copyed->max_row < M->row_dim) {
      printf("error in COPY_CRS1\n");
      printf("Nedd more Copyed->row_dim at least %d\n", M->row_dim);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "error in COPY_CRS1\n");
      fprintf(file, "Nedd more Copyed->row_dim at least %d\n", M->row_dim);
      fclose(file);
      exit(1);
   }
   if (Copyed->max_val < M->Row[M->row_dim]) {
      printf("error in COPY_CRS1\n");
      printf("Nedd more Copyed->max_val at least %ld\n", M->Row[M->row_dim]);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file,"error in COPY_CRS1\n");
      fprintf(file,"Nedd more Copyed->max_val at least %ld\n", M->Row[M->row_dim]);
      fclose(file);
      exit(1);
   }
   
#pragma omp parallel for private (j) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      for (j = M->Row[i]; j < M->Row[i+1]; j++) {
         Copyed->Col[j] = M->Col[j];
         Copyed->Val[j] = M->Val[j];
      }
      Copyed->Row[i+1] = M->Row[i+1];
   }
   
   Copyed->row_dim = M->row_dim;
   Copyed->col_dim = M->col_dim;
   
}
