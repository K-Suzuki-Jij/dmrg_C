//Ascending order

void INSERTION_SORT_LINT1_INT4_DOUBLE1(long *Int_Sort,
                                       int *Int_Array2,
                                       int *Int_Array3,
                                       int *Int_Array4,
                                       int *Int_Array5,
                                       double *D_Array2,
                                       long dim
                                       ) {
   long i,j;
   long lint_temp;
   int int_temp;
   double temp_d2;

   for (i = 1; i < dim; i++) {

      j = i;
      
      while (j > 0 && Int_Sort[j - 1] > Int_Sort[j]) {
         lint_temp     = Int_Sort[j-1];
         Int_Sort[j-1] = Int_Sort[j];
         Int_Sort[j]   = lint_temp;
         
         int_temp        = Int_Array2[j-1];
         Int_Array2[j-1] = Int_Array2[j];
         Int_Array2[j]   = int_temp;
         
         int_temp        = Int_Array3[j-1];
         Int_Array3[j-1] = Int_Array3[j];
         Int_Array3[j]   = int_temp;
         
         int_temp        = Int_Array4[j-1];
         Int_Array4[j-1] = Int_Array4[j];
         Int_Array4[j]   = int_temp;
         
         int_temp        = Int_Array5[j-1];
         Int_Array5[j-1] = Int_Array5[j];
         Int_Array5[j]   = int_temp;
         
         temp_d2       = D_Array2[j-1];
         D_Array2[j-1] = D_Array2[j];
         D_Array2[j]   = temp_d2;
         
         j--;
      }

     
   }
   
   
}
