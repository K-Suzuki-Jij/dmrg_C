//
//  EXPECTATION_VALUES.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/21.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_VALUES(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, double *Vec, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int c1 = (Dmrg_Status->sweep_now == 0);
   int c2 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 2);
   int c3 = (Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4 != Model->tot_site);
   int c4 = (Dmrg_Status->LL_site != Dmrg_Status->RR_site);

   if(c1 || c2 || c3 || c4){
      return;
   }
   
   int dim_LLLRRRRL = Dmrg_Status->dim_LLLRRRRL;
      
   double *T_V1  = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);
   double *T_V2  = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);

   double *Sz    = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *Sx    = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzSz  = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxSx  = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *Sz_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *Sx_CF = GET_ARRAY_DOUBLE1(Model->tot_site);

   OUTPUT_ENERGY(Model, Dmrg_Status);
   DMRG_EXPECTATION_ONSITE(System->Sz  , System->Sz_On  , Enviro->Sz  , Sz  , Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   DMRG_EXPECTATION_ONSITE(System->SzSz, System->SzSz_On, Enviro->SzSz, SzSz, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   DMRG_EXPECTATION_ONSITE(System->SxSx, System->SxSx_On, Enviro->SxSx, SxSx, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   
   OUTPUT_ONSITE_VALUES(Sz  , "Sz"  , Model, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SzSz, "SzSz", Model, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SxSx, "SxSx", Model, Dmrg_Status);
   
   DMRG_EXPECTATION_INTERSITE_Q0(System->Sz_CF, System->Sz, System->Sz_On, Enviro->Sz, Sz_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(Sz_CF, Sz, "Sz_CF", Model, Dmrg_Status);

   FREE_ARRAY_DOUBLE1(T_V1 );
   FREE_ARRAY_DOUBLE1(T_V2 );
   
   DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis = GET_WHOLE_BASIS_SUPERBLOCK(System, Enviro, Model, Dmrg_Basis, Dmrg_Status);
   double **T_V1_Q2  = GET_ARRAY_DOUBLE2(2, Dmrg_W_Basis->max_dim);
   double **T_V2_Q2  = GET_ARRAY_DOUBLE2(2, Dmrg_W_Basis->max_dim);
   
   EXPECTATION_INTERSITE_Q2(System->Sx_CF, System->Sx, System->Sx_On, Enviro->Sx, Sx_CF, Model->cf_origin, Vec, Model->tot_sz, T_V1_Q2, T_V2_Q2, Model->p_threads, Dmrg_W_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(Sx_CF, Sx, "Sx_CF", Model, Dmrg_Status);
   
   FREE_ARRAY_DOUBLE2(T_V1_Q2, 2);
   FREE_ARRAY_DOUBLE2(T_V2_Q2, 2);
   FREE_ARRAY_DOUBLE1(Sz  );
   FREE_ARRAY_DOUBLE1(Sx  );
   FREE_ARRAY_DOUBLE1(SzSz);
   FREE_ARRAY_DOUBLE1(SxSx);
   FREE_ARRAY_DOUBLE1(Sz_CF);
   FREE_ARRAY_DOUBLE1(Sx_CF);
   FREE_WHOLE_BASIS_SUPERBLOCK(Dmrg_W_Basis, Model);
   
}
