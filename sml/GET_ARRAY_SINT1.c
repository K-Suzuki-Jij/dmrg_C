#include <stdio.h>
#include <stdlib.h>

short int *GET_ARRAY_SINT1(long num) {
  
  long i;
  
  short int *Matrix = malloc(sizeof(short int)*num);
  
  if (Matrix == NULL) {
    printf("Error in GET_ARRAY_SINT1\n");
    printf("Need more memory(num=%ld)\n", num);
    exit(1);
  }
  
  for (i = 0; i < num; i++) {
    Matrix[i] = 0;
  }

  return Matrix;
  
}
