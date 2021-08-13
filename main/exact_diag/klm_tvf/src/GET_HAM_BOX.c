//
//  GET_HAM_BOX.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/25.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

HAM_BOX *GET_HAM_BOX(MODEL_1DKLM_TVF *Model) {
   
   int dim_onsite = Model->dim_onsite;
   int spin_loc   = Model->spin_loc;
   
   HAM_BOX *Box = malloc(sizeof(*Box));
   
   Box->Ham_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SpL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SmL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->Even_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->Even_D_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->Odd_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->Odd_D_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->Zero_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   
   
   ONSITE_MAKE_HAM_SPIN_INV_BASIS_KLM(Model, Box->Ham_On);
   ONSITE_SPL_SPIN_INV_BASIS_KLM        (spin_loc, Box->SpL_On    , 1.0);
   ONSITE_SML_SPIN_INV_BASIS_KLM        (spin_loc, Box->SmL_On    , 1.0);
   ONSITE_SZL_SPIN_INV_BASIS_KLM        (spin_loc, Box->SzL_On    , 1.0);
   ONSITE_EVEN_SPIN_INV_BASIS_KLM       (spin_loc, Box->Even_On   , 1.0);
   ONSITE_EVEN_DAGGER_SPIN_INV_BASIS_KLM(spin_loc, Box->Even_D_On , 1.0);
   ONSITE_ODD_SPIN_INV_BASIS_KLM        (spin_loc, Box->Odd_On    , 1.0);
   ONSITE_ODD_DAGGER_SPIN_INV_BASIS_KLM (spin_loc, Box->Odd_D_On  , 1.0);
   ONSITE_DIAG_SPIN_INV_BASIS_KLM       (spin_loc, Box->Zero_On   , 0.0);
   
   //For expectation values
   Box->SxL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SxC_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzC_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzLSzL_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SxLSxL_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzCSzC_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->NC_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->NCNC_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);

   ONSITE_SXL_SPIN_INV_BASIS_KLM   (spin_loc, Box->SxL_On   , 1.0);
   ONSITE_SXC_SPIN_INV_BASIS_KLM   (spin_loc, Box->SxC_On   , 1.0);
   ONSITE_SZC_SPIN_INV_BASIS_KLM   (spin_loc, Box->SzC_On   , 1.0);
   ONSITE_SZLSZL_SPIN_INV_BASIS_KLM(spin_loc, Box->SzLSzL_On, 1.0);
   ONSITE_SXLSXL_SPIN_INV_BASIS_KLM(spin_loc, Box->SxLSxL_On, 1.0);
   ONSITE_SZCSZC_SPIN_INV_BASIS_KLM(spin_loc, Box->SzCSzC_On, 1.0);
   ONSITE_NC_SPIN_INV_BASIS_KLM    (spin_loc, Box->NC_On    , 1.0);
   ONSITE_NCNC_SPIN_INV_BASIS_KLM  (spin_loc, Box->NCNC_On  , 1.0);
   
   //For SC Correlations
   Box->CCSL_On = GET_CRS2(Model->dim_ccsl_onsite, dim_onsite, 1);
   Box->CSL_On  = GET_CRS2(Model->dim_csl_onsite , dim_onsite, 2);
   
   int num;
   for (num = 0; num < Model->dim_ccsl_onsite; num++) {
      ONSITE_CCSL_KLM(num, Model->spin_loc, Box->CCSL_On[num], 1.0);
   }
   
   for (num = 0; num < Model->dim_csl_onsite; num++) {
      ONSITE_CSL_KLM(num, Model->spin_loc, Box->CSL_On[num], 1.0);
   }
   
   return Box;
   
}
