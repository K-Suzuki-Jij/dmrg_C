#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "SML.h"

void CHECK_CRS1(CRS1 *M) {
   
   long i,j;
   FILE *file;
   
   for (i = 0; i < M->row_dim; i++) {
      for (j = M->Row[i]; j < M->Row[i+1] - 1; j++){
         if (M->Col[j] >= M->Col[j+1]) {
            printf("Error in CHECK_CRS1\n");
            printf("The input Matrix is not CRS1\n");
            mkdir("SML_out",0777);
            if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
               printf("Can't open file.\n");
               exit(1);
            }
            fprintf(file,"Error in CHECK_CRS1\n");
            fprintf(file,"The input Matrix is not CRS1\n");
            fclose(file);
            exit(1);
         }
      }
   }
   
}

