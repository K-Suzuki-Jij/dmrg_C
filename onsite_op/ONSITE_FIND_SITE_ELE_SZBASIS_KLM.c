#include <stdio.h>
#include <stdlib.h>

int ONSITE_FIND_SITE_ELE_SZBASIS_KLM(int basis, int spin_loc) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
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
      printf("Error in ONSITE_FIND_SITE_ELE_SZBASIS_KLM\n");
      exit(1);
   }
   
}
