#include <stdlib.h>
#include "model.h"
#include "SML.h"

void FREE_SC_MAT_BASIS_1DHUBBARD_VF(SC_MAT_1DHUBBARD_VF *Sc_Mat, MODEL_1DHUBBARD_VF *Model) {
   
   FREE_ARRAY_CHAR2(Sc_Mat->Row_Name, Sc_Mat->max_dim);
   FREE_ARRAY_CHAR2(Sc_Mat->Col_Name, Sc_Mat->max_dim);
   FREE_ARRAY_INT1(Sc_Mat->CC_1_Num );
   FREE_ARRAY_INT1(Sc_Mat->CC_2_Num );
   FREE_ARRAY_INT1(Sc_Mat->CC_Sz    );
   FREE_ARRAY_INT1(Sc_Mat->C_Num1   );
   FREE_ARRAY_INT1(Sc_Mat->C_Num2   );
   FREE_ARRAY_INT1(Sc_Mat->C_Sz     );
   FREE_ARRAY_DOUBLE3(Sc_Mat->Mat, Model->tot_site, Sc_Mat->max_dim);
   FREE_ARRAY_DOUBLE1(Sc_Mat->F_Norm);
   free(Sc_Mat);
   
}
