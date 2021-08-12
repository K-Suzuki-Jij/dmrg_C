#include <stdio.h>
#include <stdlib.h>

char *GET_ARRAY_CHAR1(long num) {
  
  long i;
  
  char *Matrix = malloc(sizeof(char)*num);
  
  if (Matrix == NULL) {
    printf("Error in GET_ARRAY_CHAR1\n");
    printf("Need more memory(num=%ld)\n", num);
    exit(1);
  }
  
  for (i = 0; i < num; i++) {
    Matrix[i] = 0;
  }

  return Matrix;
  
}
