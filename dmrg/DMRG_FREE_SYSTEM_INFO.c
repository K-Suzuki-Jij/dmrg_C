#include <stdlib.h>
#include "dmrg.h"

void DMRG_FREE_SYSTEM_INFO(DMRG_SYSTEM_INFO *Info) {
   
   FREE_CCS1(Info->Trans_Matrix);
   FREE_CRS1(Info->Trans_Matrix_Dagger);
   FREE_ARRAY_DOUBLE1(Info->Val_DM_Dist);
   FREE_ARRAY_INT1(Info->Q_Number1);
   FREE_ARRAY_INT1(Info->Q_Number2);
   FREE_ARRAY_INT1(Info->Q_Number3);

   free(Info);
   
}
