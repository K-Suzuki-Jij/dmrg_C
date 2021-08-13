
#include <stdlib.h>
#include "model.h"
#include "SML.h"

SC_MAT_1DTKLM_VF *GET_SC_MAT_BASIS_1DTKLM_VF(MODEL_1DTKLM_VF *Model) {
   
   SC_MAT_1DTKLM_VF *Sc_Mat = malloc(sizeof(SC_MAT_1DTKLM_VF));
   Sc_Mat->max_dim          = Model->dim_ccsl_onsite;
   Sc_Mat->Row_Name         = GET_ARRAY_CHAR2(Sc_Mat->max_dim, 100);
   Sc_Mat->Col_Name         = GET_ARRAY_CHAR2(Sc_Mat->max_dim, 100);
   Sc_Mat->CCSL_Num         = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->CCSL_Sz          = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->CCSL_Ele_1       = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->CCSL_Ele_2       = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->Mat              = GET_ARRAY_DOUBLE3(Model->tot_site, Sc_Mat->max_dim, Sc_Mat->max_dim);
   
   return Sc_Mat;
   
}
