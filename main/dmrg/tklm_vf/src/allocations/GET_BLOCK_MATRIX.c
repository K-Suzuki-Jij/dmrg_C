//
//  GET_BLOCK_MATRIX.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/12.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

BLOCK *GET_BLOCK_MATRIX(MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   BLOCK *Block = malloc(sizeof(*Block));
   
   int tot_site   = Model->tot_site;
   int cf_length  = Model->tot_site/2 - 1 - Model->cf_origin;
   int max_dim    = Dmrg_Param->max_dim_system;
   int dim_onsite = Model->dim_onsite;
   int elem_num   = max_dim*max_dim*Dmrg_Param->sp_LL;
   
   if (cf_length <= 0) {
      printf("Error in GET_BLOCK_MATRIX\n");
      exit(1);
   }
   
   Block->Ham          = GET_CRS2(tot_site, max_dim, elem_num);
   Block->SpL_RE       = GET_CRS2(tot_site, max_dim, elem_num);
   Block->SmL_RE       = GET_CRS2(tot_site, max_dim, elem_num);
   Block->SzL_RE       = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CUp_1_RE     = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CDown_1_RE   = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CUp_1_D_RE   = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CDown_1_D_RE = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CUp_2_RE     = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CDown_2_RE   = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CUp_2_D_RE   = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CDown_2_D_RE = GET_CRS2(tot_site, max_dim, elem_num);
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      Block->SpL_LE       = GET_CRS2(tot_site, max_dim, elem_num);
      Block->SmL_LE       = GET_CRS2(tot_site, max_dim, elem_num);
      Block->SzL_LE       = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CUp_1_LE     = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CDown_1_LE   = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CUp_1_D_LE   = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CDown_1_D_LE = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CUp_2_LE     = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CDown_2_LE   = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CUp_2_D_LE   = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CDown_2_D_LE = GET_CRS2(tot_site, max_dim, elem_num);
   }
   
   Block->Ham_On       = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SzL_On       = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SpL_On       = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SmL_On       = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CUp_1_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CDown_1_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CUp_1_D_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CDown_1_D_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CUp_2_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CDown_2_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CUp_2_D_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CDown_2_D_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   
   Block->Tot_Sz    = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Tot_Ele   = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Tot_Ele_1 = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Tot_Ele_2 = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Dim       = GET_ARRAY_INT1(tot_site);
   
   
   ONSITE_SZL_SZBASIS_TKLM           (Model->spin_loc ,Block->SzL_On      , 1.0);
   ONSITE_SPL_SZBASIS_TKLM           (Model->spin_loc ,Block->SpL_On      , 1.0);
   ONSITE_SML_SZBASIS_TKLM           (Model->spin_loc ,Block->SmL_On      , 1.0);
   ONSITE_CUP_1_SZBASIS_TKLM         (Model->spin_loc ,Block->CUp_1_On    , 1.0);
   ONSITE_CDOWN_1_SZBASIS_TKLM       (Model->spin_loc ,Block->CDown_1_On  , 1.0);
   ONSITE_CUP_1_DAGGER_SZBASIS_TKLM  (Model->spin_loc ,Block->CUp_1_D_On  , 1.0);
   ONSITE_CDOWN_1_DAGGER_SZBASIS_TKLM(Model->spin_loc ,Block->CDown_1_D_On, 1.0);
   ONSITE_CUP_2_SZBASIS_TKLM         (Model->spin_loc ,Block->CUp_2_On    , 1.0);
   ONSITE_CDOWN_2_SZBASIS_TKLM       (Model->spin_loc ,Block->CDown_2_On  , 1.0);
   ONSITE_CUP_2_DAGGER_SZBASIS_TKLM  (Model->spin_loc ,Block->CUp_2_D_On  , 1.0);
   ONSITE_CDOWN_2_DAGGER_SZBASIS_TKLM(Model->spin_loc ,Block->CDown_2_D_On, 1.0);
   ONSITE_MAKE_HAM_SZBASIS_TKLM      (Model, Block->Ham_On);
   
   //For Expectation Values
   Block->SzC_1_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SzC_2_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SxC_1_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SxC_2_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SxL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->NC_1_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->NC_2_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SC_1SL_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SC_2SL_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);

   ONSITE_SZC_1_SZBASIS_TKLM(Model->spin_loc, Block->SzC_1_On   , 1.0);
   ONSITE_SZC_2_SZBASIS_TKLM(Model->spin_loc, Block->SzC_2_On   , 1.0);
   ONSITE_SXC_1_SZBASIS_TKLM(Model->spin_loc, Block->SxC_1_On   , 1.0);
   ONSITE_SXC_2_SZBASIS_TKLM(Model->spin_loc, Block->SxC_2_On   , 1.0);
   ONSITE_SXL_SZBASIS_TKLM(Model->spin_loc, Block->SxL_On       , 1.0);
   ONSITE_NC_1_SZBASIS_TKLM (Model->spin_loc, Block->NC_1_On    , 1.0);
   ONSITE_NC_2_SZBASIS_TKLM (Model->spin_loc, Block->NC_2_On    , 1.0);
   ONSITE_SC_1SL_SZBASIS_TKLM (Model->spin_loc, Block->SC_1SL_On, 1.0);
   ONSITE_SC_2SL_SZBASIS_TKLM (Model->spin_loc, Block->SC_2SL_On, 1.0);

   Block->TM   = malloc(sizeof(CCS1*)*tot_site/2);
   Block->TM_D = malloc(sizeof(CRS1*)*tot_site/2);

   Block->Basis_LL_LLLR  = malloc(sizeof(short*)*tot_site/2);
   Block->Basis_LR_LLLR  = malloc(sizeof(short*)*tot_site/2);
   Block->Basis_Inv_LLLR = GET_ARRAY_INT3(tot_site/2, max_dim, dim_onsite);
   Block->Dim_LLLR       = GET_ARRAY_INT1(tot_site/2);

   return Block;
}
