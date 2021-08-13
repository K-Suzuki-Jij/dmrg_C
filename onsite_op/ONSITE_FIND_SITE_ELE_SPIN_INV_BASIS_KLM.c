#include <stdio.h>
#include <stdlib.h>

int ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(int basis, int spin_loc) {
   
   ///////////////////////////////////////
   // # <->  [Cherge   ] -- (N,  P)
   // 0 <->  [         ] -- (0,  0)
   // 1 <->  [even     ] -- (1,  0)
   // 2 <->  [odd      ] -- (1,  1)
   // 3 <->  [even, odd] -- (2,  1)
   ///////////////////////////////////////
   
   int basis_charge = basis/(spin_loc + 1);
   
   if (basis_charge == 0) {
      return 0;
   }
   else if (basis_charge == 1 || basis_charge == 2) {
      return 1;
   }
   else if (basis_charge == 3) {
      return 2;
   }
   else {
      printf("Error in ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM\n");
      exit(1);
   }
   
}
