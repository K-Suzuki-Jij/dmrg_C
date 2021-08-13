//
//  EXPECTATION_SC_CORRELATIONS.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/29.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_SC_CORRELATIONS(MODEL_1DKLM_TVF *Model, EXACT_HAM_INFO *Ham_Info, EXACT_BASIS_INFO *Basis_Info, EXACT_TIME *Time) {
   
   Time->sc = omp_get_wtime();
   
   EXACT_WHOLE_BASIS_Q2 *W_Basis = GET_WHOLE_BASIS_Q2(Model, Time);
   SC_MAT_1DKLM_TVF *Sc_Mat      = GET_SC_MAT_BASIS_1DKLM_TVF(Model);

   
   SC_CORRELATIONS(0, 1, Ham_Info->Vector[0], Sc_Mat, W_Basis, Model);
   OUTPUT_SC_CORRELATIONS(Sc_Mat, 1, Model->tot_site - 1, "LSC_evev", "ScCorrelations", "EVEV", Model);
   
   SC_CORRELATIONS(1, 1, Ham_Info->Vector[0], Sc_Mat, W_Basis, Model);
   OUTPUT_SC_CORRELATIONS(Sc_Mat, 1, Model->tot_site - 1, "LSC_odod", "ScCorrelations", "ODOD", Model);
   
   FREE_SC_MAT_BASIS_1DKLM_TVF(Sc_Mat, Model);
   FREE_WHOLE_BASIS_Q2(W_Basis);
   
   Time->sc = omp_get_wtime() - Time->sc;
   
}
