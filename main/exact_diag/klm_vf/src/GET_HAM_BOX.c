//
//  GET_HAM_BOX.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/09.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

HAM_BOX *GET_HAM_BOX(MODEL_1DKLM_VF *Model) {
   
   int dim_onsite = Model->dim_onsite;
   int spin_loc   = Model->spin_loc;
   
   HAM_BOX *Box = malloc(sizeof(*Box));

   Box->Ham_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SpL_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SmL_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzL_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->CUp_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->CUp_D_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->CDown_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->CDown_D_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->Zero_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);

   ONSITE_MAKE_HAM_SZBASIS_KLM(Model, Box->Ham_On);
   ONSITE_SPL_SZBASIS_KLM         (spin_loc, Box->SpL_On    , 1.0);
   ONSITE_SML_SZBASIS_KLM         (spin_loc, Box->SmL_On    , 1.0);
   ONSITE_SZL_SZBASIS_KLM         (spin_loc, Box->SzL_On    , 1.0);
   ONSITE_CUP_SZBASIS_KLM         (spin_loc, Box->CUp_On    , 1.0);
   ONSITE_CUP_DAGGER_SZBASIS_KLM  (spin_loc, Box->CUp_D_On  , 1.0);
   ONSITE_CDOWN_SZBASIS_KLM       (spin_loc, Box->CDown_On  , 1.0);
   ONSITE_CDOWN_DAGGER_SZBASIS_KLM(spin_loc, Box->CDown_D_On, 1.0);
   ONSITE_DIAG_SZBASIS_KLM        (spin_loc, Box->Zero_On   , 0.0);

   //For expectation values
   Box->SxL_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SxC_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzC_On    = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzLSzL_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SxLSxL_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->SzCSzC_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->NC_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Box->NCNC_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   
   ONSITE_SXL_SZBASIS_KLM   (spin_loc, Box->SxL_On   , 1.0);
   ONSITE_SXC_SZBASIS_KLM   (spin_loc, Box->SxC_On   , 1.0);
   ONSITE_SZC_SZBASIS_KLM   (spin_loc, Box->SzC_On   , 1.0);
   ONSITE_SZLSZL_SZBASIS_KLM(spin_loc, Box->SzLSzL_On, 1.0);
   ONSITE_SXLSXL_SZBASIS_KLM(spin_loc, Box->SxLSxL_On, 1.0);
   ONSITE_SZCSZC_SZBASIS_KLM(spin_loc, Box->SzCSzC_On, 1.0);
   ONSITE_NC_SZBASIS_KLM    (spin_loc, Box->NC_On    , 1.0);
   ONSITE_NCNC_SZBASIS_KLM  (spin_loc, Box->NCNC_On  , 1.0);
   
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
