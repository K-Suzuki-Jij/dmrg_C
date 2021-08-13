#include <stdio.h>
#include <stdlib.h>

int ONSITE_FIND_SITE_SZ_SZBASIS_KLM(int basis, int spin_loc) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int dim_spin     = spin_loc + 1;
   int basis_loc    = basis%dim_spin;
   int basis_charge = basis/dim_spin;
   int sz_loc       = spin_loc - 2*basis_loc;
   int sz_charge;
   
   if (basis_charge == 0 || basis_charge == 3) {
      sz_charge = 0;
   }
   else if (basis_charge == 1) {
      sz_charge = 1;
   }
   else if (basis_charge == 2) {
      sz_charge = -1;
   }
   else {
      printf("Error in ONSITE_FIND_SITE_SZ_SZBASIS_KLM\n");
      exit(1);
   }
   
   return sz_charge + sz_loc;
   
}
