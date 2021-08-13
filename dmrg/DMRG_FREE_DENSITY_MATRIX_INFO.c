#include <stdlib.h>
#include "dmrg.h"

void DMRG_FREE_DENSITY_MATRIX_INFO(DMRG_DENSITY_MATRIX_INFO *Info) {
   
   FREE_ARRAY_INT1(Info->Dim_Block);
   FREE_ARRAY_INT1(Info->Q_Number1_Block);
   FREE_ARRAY_INT1(Info->Q_Number2_Block);
   FREE_ARRAY_INT1(Info->Q_Number3_Block);
   free(Info);

}
