#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SML.h"

void WRITE_CRS1(CRS1 *M, char Name[]) {
  
  FILE *file;
  char s1[] = "_Val";
  char s2[] = "_Col";
  char s3[] = "_Row";
  char s4[] = "_row_dim";
  char s5[] = "_col_dim";
  char s6[] = "_val_num";

  char Name1[100];
  char Name2[100];
  char Name3[100];
  char Name4[100];
  char Name5[100];
  char Name6[100];
  
  strcpy(Name1, Name);
  strcpy(Name2, Name);
  strcpy(Name3, Name);
  strcpy(Name4, Name);
  strcpy(Name5, Name);
  strcpy(Name6, Name);

  strcat(Name1, s1);
  strcat(Name2, s2);
  strcat(Name3, s3);
  strcat(Name4, s4);
  strcat(Name5, s5);
  strcat(Name6, s6);

  file = fopen(Name1, "wb");
  fwrite(M->Val, sizeof(double), M->Row[M->row_dim], file);
  fclose(file);
  
  file = fopen(Name2, "wb");
  fwrite(M->Col, sizeof(int), M->Row[M->row_dim], file);
  fclose(file);
  
  file = fopen(Name3, "wb");
  fwrite(M->Row, sizeof(long), M->row_dim + 1, file);
  fclose(file);
  
  file = fopen(Name4, "wb");
  fwrite(&M->row_dim, sizeof(long), 1, file);
  fclose(file);
  
  file = fopen(Name5, "wb");
  fwrite(&M->col_dim, sizeof(long), 1, file);
  fclose(file);
  
  file = fopen(Name6, "wb");
  fwrite(&M->Row[M->row_dim], sizeof(long), 1, file);
  fclose(file);
  
}
