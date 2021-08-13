//
//  EXPECTATION_SC_CORRELATIONS.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/12.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_SC_CORRELATIONS(MODEL_1DKLM_VF *Model, EXACT_HAM_INFO *Ham_Info, EXACT_BASIS_INFO *Basis_Info, EXACT_TIME *Time) {
   
   Time->sc = omp_get_wtime();
   
   EXACT_WHOLE_BASIS_Q3 *W_Basis = GET_WHOLE_BASIS_Q3(Model, Time);
   SC_MAT_1DKLM_VF      *Sc_Mat  = GET_SC_MAT_BASIS_1DKLM_VF(Model);
   
   int del_sz;
   int max_del_sz = 4*Model->spin_loc + 2;
   char Name1[100],Name2[100];
   
   for (del_sz = -max_del_sz; del_sz <= max_del_sz; del_sz++) {
      SC_CORRELATIONS(del_sz, 1, Ham_Info->Vector[0], Sc_Mat, W_Basis, Model);
      sprintf(Name1, "LSC_%d_%d", del_sz, del_sz);
      sprintf(Name2, "%d_%d", del_sz, del_sz);
      OUTPUT_SC_CORRELATIONS(Sc_Mat, 1, Model->tot_site - 1, Name1, "ScCorrelations", Name2, Model);
   }
   
   FREE_SC_MAT_BASIS_1DKLM_VF(Sc_Mat, Model);
   FREE_WHOLE_BASIS_Q3(W_Basis, Model->tot_site*Model->spin_loc + Model->tot_ele);
   
   Time->sc = omp_get_wtime() - Time->sc;
   
}
