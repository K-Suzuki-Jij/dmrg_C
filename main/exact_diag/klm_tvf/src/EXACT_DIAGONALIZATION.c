//
//  EXACT_DIAGONALIZATION.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/25.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXACT_DIAGONALIZATION(MODEL_1DKLM_TVF *Model, EXACT_PARAMETER *Param, EXACT_BASIS_INFO *Basis_Info, EXACT_HAM_INFO *Ham_Info, EXACT_TIME *Time) {
   
   Time->total = omp_get_wtime();
   
   Ham_Info->Ham = GET_HAM(Model, Basis_Info, Param, Time);
   EXACT_DIAGONALIZE_HAMILTONIAN(Ham_Info, Param, Time, Model->p_threads);
   FREE_CRS1(Ham_Info->Ham);
   
   OUTPUT_ENERGY(Model, Ham_Info);
   
   EXPECTATION_VALUES(Model, Ham_Info, Basis_Info, Time);
   
   EXPECTATION_SC_CORRELATIONS(Model, Ham_Info, Basis_Info, Time);
   
   Time->total = omp_get_wtime() - Time->total;
   
   PRINT_STATUS(Model, Time, Ham_Info, Param);
   
}
