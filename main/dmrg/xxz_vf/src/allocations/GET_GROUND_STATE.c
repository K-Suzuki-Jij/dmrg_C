//
//  GET_GROUND_STATE.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/27.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void GET_GROUND_STATE(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Status->GS_Vec = GET_ARRAY_DOUBLE1(Dmrg_Basis->dim_LLLRRRRL);
   
   if (Dmrg_Basis->dim_LLLRRRRL < Dmrg_Param->dim_LLLRRRRL_limit) {
      CRS1 *Ham = GET_HAM_LLLRRRRL(Dmrg_Basis, System, Enviro, Model, Dmrg_Status, Dmrg_Time);
      DMRG_DIAGONALIZE_SUPERBLOCK_SYM(Ham, Dmrg_Time, Dmrg_Param, Dmrg_Status, Model->p_threads);
      FREE_CRS1(Ham);
   }
   else if (Dmrg_Basis->dim_LLLRRRRL >= Dmrg_Param->dim_LLLRRRRL_limit && (strcmp(Model->BC, "OBC") == 0 || strcmp(Model->BC, "SSD") == 0)) {
      DIAGONALIZE_SUPERBLOCK(Dmrg_Basis, System, Enviro, Model, Dmrg_Param, Dmrg_Status, Dmrg_Time);
   }
   else {
      printf("Error in RENORMALIZE\n");
      printf("dim_LLLRRRRL = %d >= %d\n", Dmrg_Basis->dim_LLLRRRRL, Dmrg_Param->dim_LLLRRRRL_limit);
      printf("BC = %s\n", Model->BC);
      exit(1);
   }

   
}
