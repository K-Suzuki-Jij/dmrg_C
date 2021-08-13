#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double DMRG_SSD_COEFF(int LL_site, int tot_site, char BC[], char Block_Name[], char Inter_Name[]) {
   
   int site;
   
   if(strcmp(Block_Name,"LL") == 0){
      site = LL_site + 1;
   }
   else if(strcmp(Block_Name,"LR") == 0){
      site = LL_site + 2;
   }
   else if(strcmp(Block_Name,"RL") == 0){
      site = LL_site + 3; //tot_site - (RR_site + 1);
   }
   else {
      printf("Error in SSD_COEFF\n");
      printf("Name=%s\n",Block_Name);
      exit(1);
   }
   
   double temp;
   
   if (strcmp(BC, "SSD") == 0) {
      if (strcmp(Inter_Name, "Onsite") == 0) {
         temp = sin(M_PI*(site - 0.5)/tot_site);
      }
      else if (strcmp(Inter_Name, "Intersite") == 0) {
         temp = sin(M_PI*site/tot_site);
      }
      else {
         printf("Error in SSD_COEFF\n");
         exit(1);
      }
   }
   else {
      temp = 1.0;
   }
   
   
   return temp*temp;
   
}
