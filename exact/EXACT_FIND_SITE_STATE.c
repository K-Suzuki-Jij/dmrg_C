#include <stdio.h>

int EXACT_FIND_SITE_STATE(long basis, int site, int dim_onsite) {
   
   int i;
   
   for (i = 0; i < site; i++) {
      basis = basis/(long)dim_onsite;
   }
   
   basis = basis%dim_onsite;
   
   return (int)basis;
   
}
