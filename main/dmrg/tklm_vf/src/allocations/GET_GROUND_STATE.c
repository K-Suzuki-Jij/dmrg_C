//
//  GET_GROUND_STATE.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/28.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void GET_GROUND_STATE(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Status->GS_Vec = GET_ARRAY_DOUBLE1(Dmrg_Basis->dim_LLLRRRRL);
   
   if (Dmrg_Basis->dim_LLLRRRRL < Dmrg_Param->dim_LLLRRRRL_limit) {
      CRS1 *Ham = GET_HAM_LLLRRRRL(Dmrg_Basis, System, Enviro, Model, Dmrg_Status, Dmrg_Time);
      FREE_ARRAY_INT1(Dmrg_Basis->Inv_LLLRRRRL);
      DMRG_DIAGONALIZE_SUPERBLOCK_SYM(Ham, Dmrg_Time, Dmrg_Param, Dmrg_Status, Model->p_threads);
      DMRG_RE_ALLOCATE_INV_LLLRRRRL(Dmrg_Basis, System->Dim[Dmrg_Status->LL_site], Enviro->Dim[Dmrg_Status->RR_site], Model->dim_onsite, Model->p_threads);
      FREE_CRS1(Ham);
   }
   else {
      printf("Error in RENORMALIZE\n");
      printf("dim_LLLRRRRL = %d >= %d\n", Dmrg_Basis->dim_LLLRRRRL, Dmrg_Param->dim_LLLRRRRL_limit);
      printf("BC = %s\n", Model->BC);
      exit(1);
   }

}
