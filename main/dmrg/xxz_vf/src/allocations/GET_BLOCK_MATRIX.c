//
//  GET_BLOCK_MATRIX.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/12.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

BLOCK *GET_BLOCK_MATRIX(MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
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
   
   Block->Ham   = GET_CRS2(tot_site, max_dim, elem_num);
   Block->Sz_RE = GET_CRS2(tot_site, max_dim, elem_num);
   Block->Sp_RE = GET_CRS2(tot_site, max_dim, elem_num);
   Block->Sm_RE = GET_CRS2(tot_site, max_dim, elem_num);
   
   if (strcmp(Model->BC, "PBC") == 0) {
      Block->Sz_LE = GET_CRS2(tot_site, max_dim, elem_num);
      Block->Sp_LE = GET_CRS2(tot_site, max_dim, elem_num);
      Block->Sm_LE = GET_CRS2(tot_site, max_dim, elem_num);
   }
   
   Block->Sz_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Sp_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Sm_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Ham_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Sx_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SzSz_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SxSx_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   
   Block->Sz    = GET_CRS2(tot_site/2, max_dim, elem_num);
   Block->SzSz  = GET_CRS2(tot_site/2, max_dim, elem_num);
   Block->Sx    = GET_CRS2(tot_site/2, max_dim, elem_num);
   Block->SxSx  = GET_CRS2(tot_site/2, max_dim, elem_num);
   Block->Sz_CF = GET_CRS2(cf_length, max_dim, elem_num);
   Block->Sx_CF = GET_CRS2(cf_length, max_dim, elem_num);
   
   Block->Tot_Sz = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Dim    = GET_ARRAY_INT1(tot_site);
   
   ONSITE_SZ_SZBASIS_HB  (Model->spin ,Block->Sz_On  , 1.0);
   ONSITE_SP_SZBASIS_HB  (Model->spin ,Block->Sp_On  , 1.0);
   ONSITE_SM_SZBASIS_HB  (Model->spin ,Block->Sm_On  , 1.0);
   ONSITE_SX_SZBASIS_HB  (Model->spin ,Block->Sx_On  , 1.0);
   ONSITE_SZSZ_SZBASIS_HB(Model->spin ,Block->SzSz_On, 1.0);
   ONSITE_SXSX_SZBASIS_HB(Model->spin ,Block->SxSx_On, 1.0);
   MAKE_ONSITE_HAM       (Model       ,Block->Ham_On );
   
   
   return Block;
}
