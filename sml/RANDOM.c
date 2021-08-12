#include <stdlib.h>

long RANDOM(long min, long max) {
   
   return min + (long)(rand()*(max - min + 1)/(1 + (long)RAND_MAX));
   
}
