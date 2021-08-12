//Descending order

long BINARY_SEARCH_LINT1(long *Array, long imin, long imax, long target_val) {
   
   if (imax < imin) {
      return -1;
   }
   else {
      long imid = imin + (imax - imin)/2;
      if (Array[imid] > target_val) {
         return BINARY_SEARCH_LINT1(Array, imin, imid - 1, target_val);
      }
      else if (Array[imid] < target_val) {
         return BINARY_SEARCH_LINT1(Array, imid + 1, imax, target_val);
      }
      else {
         return imid;
      }
   }

}
