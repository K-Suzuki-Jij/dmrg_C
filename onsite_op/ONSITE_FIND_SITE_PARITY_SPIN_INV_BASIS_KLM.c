#include <stdio.h>
#include <stdlib.h>

int ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(int basis, int spin_loc) {
   
   ///////////////////////////////////////
   // # <->  [Cherge   ] -- (N,  P)
   // 0 <->  [         ] -- (0,  0)
   // 1 <->  [even     ] -- (1,  0)
   // 2 <->  [odd      ] -- (1,  1)
   // 3 <->  [even, odd] -- (2,  1)
   ///////////////////////////////////////
   
   int dim_spin     = spin_loc + 1;
   int basis_loc    = basis%dim_spin;
   int basis_charge = basis/dim_spin;
   int dim_parity   = 2;
   int parity_loc,parity_charge;
   
   if (basis_charge == 0 || basis_charge == 1) {
      parity_charge = 0;
   }
   else if (basis_charge == 2 || basis_charge == 3) {
      parity_charge = 1;
   }
   else {
      printf("Error in ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM\n");
      exit(1);
   }
   
   parity_loc = basis_loc%dim_parity;
   
   return (parity_loc + parity_charge)%2;
   
}
