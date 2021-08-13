//
//  GET_BLOCK_MATRIX.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

BLOCK *GET_BLOCK_MATRIX(MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   BLOCK *Block = malloc(sizeof(*Block));
   
   int tot_site   = Model->tot_site;
   int max_dim    = Dmrg_Param->max_dim_system;
   int dim_onsite = Model->dim_onsite;
   int elem_num   = max_dim*max_dim*Dmrg_Param->sp_LL;

   Block->Ham       = GET_CRS2(tot_site, max_dim, elem_num);
   Block->SpL_RE    = GET_CRS2(tot_site, max_dim, elem_num);
   Block->SmL_RE    = GET_CRS2(tot_site, max_dim, elem_num);
   Block->SzL_RE    = GET_CRS2(tot_site, max_dim, elem_num);
   Block->Even_RE   = GET_CRS2(tot_site, max_dim, elem_num);
   Block->Odd_RE    = GET_CRS2(tot_site, max_dim, elem_num);
   Block->Even_D_RE = GET_CRS2(tot_site, max_dim, elem_num);
   Block->Odd_D_RE  = GET_CRS2(tot_site, max_dim, elem_num);
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      Block->SpL_LE    = GET_CRS2(tot_site, max_dim, elem_num);
      Block->SmL_LE    = GET_CRS2(tot_site, max_dim, elem_num);
      Block->SzL_LE    = GET_CRS2(tot_site, max_dim, elem_num);
      Block->Even_LE   = GET_CRS2(tot_site, max_dim, elem_num);
      Block->Odd_LE    = GET_CRS2(tot_site, max_dim, elem_num);
      Block->Even_D_LE = GET_CRS2(tot_site, max_dim, elem_num);
      Block->Odd_D_LE  = GET_CRS2(tot_site, max_dim, elem_num);
   }
   
   Block->Tot_Parity = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Tot_Ele    = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Dim        = GET_ARRAY_INT1(tot_site);
   
   Block->SzL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SpL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SmL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Even_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Odd_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Even_D_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Odd_D_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Ham_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   
   ONSITE_SZL_SPIN_INV_BASIS_KLM        (Model->spin_loc ,Block->SzL_On   , 1.0);
   ONSITE_SPL_SPIN_INV_BASIS_KLM        (Model->spin_loc ,Block->SpL_On   , 1.0);
   ONSITE_SML_SPIN_INV_BASIS_KLM        (Model->spin_loc ,Block->SmL_On   , 1.0);
   ONSITE_EVEN_SPIN_INV_BASIS_KLM       (Model->spin_loc ,Block->Even_On  , 1.0);
   ONSITE_ODD_SPIN_INV_BASIS_KLM        (Model->spin_loc ,Block->Odd_On   , 1.0);
   ONSITE_EVEN_DAGGER_SPIN_INV_BASIS_KLM(Model->spin_loc ,Block->Even_D_On, 1.0);
   ONSITE_ODD_DAGGER_SPIN_INV_BASIS_KLM (Model->spin_loc ,Block->Odd_D_On , 1.0);
   ONSITE_MAKE_HAM_SPIN_INV_BASIS_KLM   (Model           ,Block->Ham_On);
   
   //For Expectation Values
   Block->SxL_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SzC_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SxC_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SCSL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->NC_On      = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   
   ONSITE_SXL_SPIN_INV_BASIS_KLM (Model->spin_loc, Block->SxL_On , 1.0);
   ONSITE_SZC_SPIN_INV_BASIS_KLM (Model->spin_loc, Block->SzC_On , 1.0);
   ONSITE_SXC_SPIN_INV_BASIS_KLM (Model->spin_loc, Block->SxC_On , 1.0);
   ONSITE_SCSL_SPIN_INV_BASIS_KLM(Model->spin_loc, Block->SCSL_On, 1.0);
   ONSITE_NC_SPIN_INV_BASIS_KLM  (Model->spin_loc, Block->NC_On  , 1.0);
   
   Block->TM   = malloc(sizeof(CCS1*)*tot_site/2);
   Block->TM_D = malloc(sizeof(CRS1*)*tot_site/2);
   
   Block->Basis_LL_LLLR  = malloc(sizeof(short*)*tot_site/2);
   Block->Basis_LR_LLLR  = malloc(sizeof(short*)*tot_site/2);
   Block->Basis_Inv_LLLR = GET_ARRAY_INT3(tot_site/2, max_dim, dim_onsite);
   Block->Dim_LLLR       = GET_ARRAY_INT1(tot_site/2);
   
   return Block;
   
}
