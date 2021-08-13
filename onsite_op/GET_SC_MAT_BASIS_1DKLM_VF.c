
#include <stdlib.h>
#include "model.h"
#include "SML.h"

SC_MAT_1DKLM_VF *GET_SC_MAT_BASIS_1DKLM_VF(MODEL_1DKLM_VF *Model) {
   
   SC_MAT_1DKLM_VF *Sc_Mat = malloc(sizeof(SC_MAT_1DKLM_VF));
   Sc_Mat->max_dim  = Model->dim_lspin*Model->dim_lspin*2*Model->dim_lspin*Model->dim_lspin*2 + 2*Model->dim_lspin*Model->dim_lspin;
   Sc_Mat->Row_Name = GET_ARRAY_CHAR2(Sc_Mat->max_dim, 100);
   Sc_Mat->Col_Name = GET_ARRAY_CHAR2(Sc_Mat->max_dim, 100);
   Sc_Mat->CC_1_Num = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->CC_2_Num = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->CC_Sz    = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->C_Num1   = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->C_Num2   = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->C_Sz     = GET_ARRAY_INT1(Sc_Mat->max_dim);
   Sc_Mat->Mat      = GET_ARRAY_DOUBLE3(Model->tot_site, Sc_Mat->max_dim, Sc_Mat->max_dim);
   Sc_Mat->F_Norm   = GET_ARRAY_DOUBLE1(Model->tot_site);
   
   return Sc_Mat;
   
}
