#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/stat.h>
#include "SML.h"

void DIAG_ADD_CRS1(CRS1 *M, double add, int p_threads) {
   
   long i,j,count;
   FILE *file;
   
   count = 0;
#pragma omp parallel for private (j) reduction (+:count) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      for (j = M->Row[i]; j < M->Row[i+1]; j++) {
         if (i == M->Col[j]) {
            M->Val[j] = M->Val[j] + add;
            count++;
            break;
         }
      }
   }
   
   if (count != M->row_dim) {
      printf("Error in DIAG_ADD_CRS1\n");
      printf("We can not add a value to the diagonal elements\n");
      mkdir("SML_out",0777);
      if((file = fopen("./SML_out/error.txt","a+")) == NULL){
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file,"Error in DIAG_ADD_CRS1\n");
      fprintf(file,"We can not add a value to the diagonal elements\n");
      fclose(file);
      exit(1);
   }
   
}
