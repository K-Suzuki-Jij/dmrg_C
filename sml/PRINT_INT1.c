#include <stdio.h>

void PRINT_INT1(int *Vec, long dim, char Name[]) {
   
   long i;
   
   for (i = 0; i < dim; i++) {
      printf("%s[%ld]=%d\n", Name, i, Vec[i]);
   }
   
}
