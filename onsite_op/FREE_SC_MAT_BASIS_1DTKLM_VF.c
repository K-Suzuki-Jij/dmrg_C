#include <stdlib.h>
#include "model.h"
#include "SML.h"

void FREE_SC_MAT_BASIS_1DTKLM_VF(SC_MAT_1DTKLM_VF *Sc_Mat, MODEL_1DTKLM_VF *Model) {
   
   FREE_ARRAY_CHAR2(Sc_Mat->Row_Name, Sc_Mat->max_dim);
   FREE_ARRAY_CHAR2(Sc_Mat->Col_Name, Sc_Mat->max_dim);
   FREE_ARRAY_INT1(Sc_Mat->CCSL_Num  );
   FREE_ARRAY_INT1(Sc_Mat->CCSL_Sz   );
   FREE_ARRAY_INT1(Sc_Mat->CCSL_Ele_1);
   FREE_ARRAY_INT1(Sc_Mat->CCSL_Ele_2);
   FREE_ARRAY_DOUBLE3(Sc_Mat->Mat, Model->tot_site, Sc_Mat->max_dim);
   free(Sc_Mat);
   
}
