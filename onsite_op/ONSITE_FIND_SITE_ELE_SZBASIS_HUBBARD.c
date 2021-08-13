#include <stdio.h>
#include <stdlib.h>

int ONSITE_FIND_SITE_ELE_SZBASIS_HUBBARD(int basis) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (basis == 0) {
      return 0;
   }
   else if (basis == 1 || basis == 2) {
      return 1;
   }
   else if (basis == 3) {
      return 2;
   }
   else {
      printf("Error in ONSITE_FIND_SITE_ELE_SZBASIS_HUBBARD\n");
      printf("basis = %d\n", basis);
      exit(1);
   }
      
}
