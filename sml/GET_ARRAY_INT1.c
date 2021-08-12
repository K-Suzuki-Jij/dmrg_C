#include <stdio.h>
#include <stdlib.h>

int *GET_ARRAY_INT1(long num) {
  
  long i;
  
  int *Matrix = malloc(sizeof(int)*num);
  
  if (Matrix == NULL) {
    printf("Error in GET_ARRAY_INT1\n");
    printf("Need more memory(num=%ld)\n", num);
    exit(1);
  }
  
  for (i = 0; i < num; i++) {
    Matrix[i] = 0;
  }

  return Matrix;
  
}
