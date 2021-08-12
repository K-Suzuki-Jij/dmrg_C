#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SML.h"

void READ_CRS1(CRS1 *M, char Name[]) {
   
   FILE *file;
   char s0[] = "_val_num";
   char s1[] = "_Val";
   char s2[] = "_Col";
   char s3[] = "_Row";
   char s4[] = "_row_dim";
   char s5[] = "_col_dim";
   
   char Name0[100];
   char Name1[100];
   char Name2[100];
   char Name3[100];
   char Name4[100];
   char Name5[100];

   strcpy(Name0, Name);
   strcpy(Name1, Name);
   strcpy(Name2, Name);
   strcpy(Name3, Name);
   strcpy(Name4, Name);
   strcpy(Name5, Name);

   strcat(Name0, s0);
   strcat(Name1, s1);
   strcat(Name2, s2);
   strcat(Name3, s3);
   strcat(Name4, s4);
   strcat(Name5, s5);
   
   long val_num;
   
   file = fopen(Name0,"rb");
   if (file == NULL) {
      printf("error int READ_CRS1\n");
      printf("cant open the file named by %s\n", Name0);
      exit(1);
   }
   fread(&val_num,sizeof(long), 1, file);
   fclose(file);
   
   file = fopen(Name1,"rb");
   if (file == NULL) {
      printf("error int READ_CRS1\n");
      printf("cant open the file named by %s\n", Name1);
      exit(1);
   }
   fread(M->Val, sizeof(double), val_num, file);
   fclose(file);
   
   file = fopen(Name2,"rb");
   if (file == NULL) {
      printf("error int READ_CRS1\n");
      printf("cant open the file named by %s\n", Name2);
      exit(1);
   }
   fread(M->Col, sizeof(int), val_num, file);
   fclose(file);
   
   file = fopen(Name3,"rb");
   if (file == NULL) {
      printf("error int READ_CRS1\n");
      printf("cant open the file named by %s\n", Name3);
      exit(1);
   }
   fread(M->Row, sizeof(long), M->row_dim + 1, file);
   fclose(file);
   
   file = fopen(Name4,"rb");
   if (file == NULL) {
      printf("error int READ_CRS1\n");
      printf("cant open the file named by %s\n", Name4);
      exit(1);
   }
   fread(&M->row_dim, sizeof(long), 1, file);
   fclose(file);
   
   file = fopen(Name5,"rb");
   if (file == NULL) {
      printf("error int READ_CRS1\n");
      printf("cant open the file named by %s\n", Name5);
      exit(1);
   }
   fread(&M->col_dim, sizeof(long), 1, file);
   fclose(file);
   
}
