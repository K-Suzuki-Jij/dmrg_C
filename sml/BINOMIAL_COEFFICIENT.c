#include <stdio.h>
#include <stdlib.h>

long BINOMIAL_COEFFICIENT(int n, int k) {
   
   if (n <= 0 || k < 0 || k > n) {
      printf("Error in BINOMIAL_COEFFICIENT\n");
      printf("n=%d,k=%d\n",n,k);
      exit(1);
   }
   
   if (k == 0 || k == n) {
      return 1;
   }
   else {
      return BINOMIAL_COEFFICIENT(n - 1, k) + BINOMIAL_COEFFICIENT(n - 1, k - 1);
   }
   
}
