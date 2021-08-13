//
//  RENORMALIZE.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void RENORMALIZE(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->total = omp_get_wtime();
   Dmrg_Status->tot_iter_now++;
   
   if (Dmrg_Status->LL_site == 0) {
      MAKE_OP_EDGE(System, Model);
   }
    
   DMRG_BASIS *Dmrg_Basis = GET_BASIS_SUPERBLOCK(System, Enviro, Model, Dmrg_Status, Dmrg_Time);

   GET_GROUND_STATE(Dmrg_Basis, System, Enviro, Model, Dmrg_Param, Dmrg_Status, Dmrg_Time);
   
   EXPECTATION_VALUES(System, Enviro, Model, Dmrg_Status->GS_Vec, Dmrg_Time, Dmrg_Basis, Dmrg_Param, Dmrg_Status);
   
   EXPECTATION_SC_CORRELATIONS(System, Enviro, Model, Dmrg_Param, Dmrg_Status->GS_Vec, Dmrg_Time, Dmrg_Basis, Dmrg_Status);
         
   DMRG_SYSTEM_INFO *Dmrg_System = DMRG_GET_SYSTEM_INFO_Q2(Dmrg_Basis->Tot_Ele_LLLR, Dmrg_Basis->Tot_Sz_LLLR, Dmrg_Basis, Dmrg_Status->GS_Vec, Dmrg_Param->max_dim_system, Model->p_threads, Dmrg_Time);
   
   OUTPUT_T_ERROR(Dmrg_System, Dmrg_Status, Model);
   
   TRANSFORM_MATRIX(System, Enviro, Dmrg_Basis, Dmrg_System, Dmrg_Status, Model, Dmrg_Time);
      
   STORE_T_MATRIX(Dmrg_Basis, Dmrg_System, System, Model, Dmrg_Status, Dmrg_Status->LL_site, Model->dim_onsite);
   
   DMRG_FREE_MEMORY(Dmrg_Basis, Dmrg_System);
   
   FREE_ARRAY_DOUBLE1(Dmrg_Status->GS_Vec);

   PRINT_STATUS(Dmrg_Status, Dmrg_Time);
   
}
