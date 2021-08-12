#include <stdio.h>
#include <stdlib.h>

long *GET_ARRAY_LINT1(long num) {
  
  long i;
  
  long *Matrix = malloc(sizeof(long)*num);
  if (Matrix == NULL) {
    printf("Error in GET_ARRAY_LINT1\n");
    printf("Need more memory(num=%ld)\n", num);
    exit(1);
  }
  
  for (i = 0; i < num; i++) {
    Matrix[i] = 0;
  }

  return Matrix;
  
}
