#include <stdio.h>

void PRINT_LINT1(long *Vec, long dim, char Name[]) {
   
   long i;
   
   for (i = 0; i < dim; i++) {
      printf("%s[%ld]=%ld\n", Name, i, Vec[i]);
   }
   
}
