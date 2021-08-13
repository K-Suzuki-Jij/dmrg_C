#include <stdio.h>
#include <stdlib.h>

int ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(int basis, int spin_loc) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int dim_spin       = spin_loc + 1;
   int dim_charge     = 4;
   int basis_loc      = basis%dim_spin;
   int basis_charge_1 = basis/(spin_loc + 1)/dim_charge;
   int basis_charge_2 = (basis/(spin_loc + 1))%dim_charge;
   int sz_loc         = spin_loc - 2*basis_loc;
   int sz_charge_1,sz_charge_2;
   
   if (basis_charge_1 == 0 || basis_charge_1 == 3) {
      sz_charge_1 = 0;
   }
   else if (basis_charge_1 == 1) {
      sz_charge_1 = 1;
   }
   else if (basis_charge_1 == 2) {
      sz_charge_1 = -1;
   }
   else {
      printf("Error in ONSITE_FIND_SITE_SZ_SZBASIS_TKLM\n");
      exit(1);
   }
   if (basis_charge_2 == 0 || basis_charge_2 == 3) {
      sz_charge_2 = 0;
   }
   else if (basis_charge_2 == 1) {
      sz_charge_2 = 1;
   }
   else if (basis_charge_2 == 2) {
      sz_charge_2 = -1;
   }
   else {
      printf("Error in ONSITE_FIND_SITE_SZ_SZBASIS_TKLM\n");
      exit(1);
   }
   
   return sz_charge_1 + sz_charge_2 + sz_loc;
   
}
