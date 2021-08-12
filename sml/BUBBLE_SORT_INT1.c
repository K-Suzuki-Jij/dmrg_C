//Descending order

void BUBBLE_SORT_INT1(int *Array, int dim) {
   
   int i,j,count,temp;

   for(i = 0; i < dim - 1; i++) {
      count = 0;
      for (j = 0; j < dim - 1 - i; j++) {
         if (Array[j] < Array[j+1]) {
            temp       = Array[j];
            Array[j]   = Array[j+1];
            Array[j+1] = temp;
            count = 1;
         }
      }
      
      if (count == 0) {
         break;
      }
   }
}
