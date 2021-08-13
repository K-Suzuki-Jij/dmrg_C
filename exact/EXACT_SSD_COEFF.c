#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


double EXACT_SSD_COEFF(int site, int tot_site, char BC[], char Inter_Name[]) {
   
   site = site + 1; //site = 0 means the left edge
   
   double temp;
   
   if (strcmp(BC, "SSD") == 0) {
      if (strcmp(Inter_Name, "Onsite") == 0) {
         temp = sin(M_PI*(site - 0.5)/tot_site);
      }
      else if (strcmp(Inter_Name, "Intersite") == 0) {
         temp = sin(M_PI*site/tot_site);
      }
      else {
         printf("Error in EXACT_SSD_COEFF\n");
         exit(1);
      }
   }
   else {
      temp = 1.0;
   }
   
   return temp*temp;

}
