//
//  MAKE_ELEMENT_HAM.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/25.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM(long basis, EXACT_A_BASIS *A_Basis, HAM_BOX *Ham_Box, MODEL_1DKLM_TVF *Model) {
   
   int site,i,sign,n_ele;
   int tot_site   = Model->tot_site;
   int dim_onsite = Model->dim_onsite;
   long elem_num  = 0;
   double ssd_coeef;
   
   
   //ONSITE
   for(site = 0; site < tot_site; site++){
      ssd_coeef = EXACT_SSD_COEFF(site, tot_site, Model->BC, "Onsite");
      EXACT_MAKE_ELEM_ON(basis, site, dim_onsite, Ham_Box->Ham_On, &elem_num, ssd_coeef, A_Basis);
   }
   
   //INTERSITE
   for(site = 0; site < tot_site - 1; site++){
      ssd_coeef = EXACT_SSD_COEFF(site, tot_site, Model->BC, "Intersite");
      EXACT_MAKE_ELEM_INTER(basis, site, site + 1, dim_onsite, Ham_Box->SzL_On, Ham_Box->SzL_On, &elem_num, Model->I_z*ssd_coeef     , 1, A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, site, site + 1, dim_onsite, Ham_Box->SpL_On, Ham_Box->SmL_On, &elem_num, 0.5*Model->I_xy*ssd_coeef, 1, A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, site, site + 1, dim_onsite, Ham_Box->SmL_On, Ham_Box->SpL_On, &elem_num, 0.5*Model->I_xy*ssd_coeef, 1, A_Basis);
      n_ele = ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(EXACT_FIND_SITE_STATE(basis, site, dim_onsite), Model->spin_loc);
      if (n_ele%2 == 1) {
         sign = 1;
      }
      else {
         sign = -1;
      }
      EXACT_MAKE_ELEM_INTER(basis, site, site + 1, dim_onsite, Ham_Box->Even_D_On, Ham_Box->Even_On  , &elem_num, Model->t*ssd_coeef, sign   , A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, site, site + 1, dim_onsite, Ham_Box->Even_On  , Ham_Box->Even_D_On, &elem_num, Model->t*ssd_coeef, -1*sign, A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, site, site + 1, dim_onsite, Ham_Box->Odd_D_On , Ham_Box->Odd_On   , &elem_num, Model->t*ssd_coeef, sign   , A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, site, site + 1, dim_onsite, Ham_Box->Odd_On   , Ham_Box->Odd_D_On , &elem_num, Model->t*ssd_coeef, -1*sign, A_Basis);
   }
   
   if (strcmp(Model->BC, "PBC") == 0) {
      EXACT_MAKE_ELEM_INTER(basis, 0, Model->tot_site - 1, dim_onsite, Ham_Box->SzL_On, Ham_Box->SzL_On, &elem_num, Model->I_z     , 1, A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, 0, Model->tot_site - 1, dim_onsite, Ham_Box->SpL_On, Ham_Box->SmL_On, &elem_num, 0.5*Model->I_xy, 1, A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, 0, Model->tot_site - 1, dim_onsite, Ham_Box->SmL_On, Ham_Box->SpL_On, &elem_num, 0.5*Model->I_xy, 1, A_Basis);
      
      n_ele = 0;
      for(site = 0; site < tot_site - 1; site++){
         n_ele = n_ele + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(EXACT_FIND_SITE_STATE(basis, site, dim_onsite), Model->spin_loc);
      }
      if (n_ele%2 == 1) {
         sign = 1;
      }
      else {
         sign = -1;
      }
      
      EXACT_MAKE_ELEM_INTER(basis, 0, Model->tot_site - 1, dim_onsite, Ham_Box->Even_D_On, Ham_Box->Even_On  , &elem_num, Model->t, sign   , A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, 0, Model->tot_site - 1, dim_onsite, Ham_Box->Even_On  , Ham_Box->Even_D_On, &elem_num, Model->t, -1*sign, A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, 0, Model->tot_site - 1, dim_onsite, Ham_Box->Odd_D_On , Ham_Box->Odd_On   , &elem_num, Model->t, sign   , A_Basis);
      EXACT_MAKE_ELEM_INTER(basis, 0, Model->tot_site - 1, dim_onsite, Ham_Box->Odd_On   , Ham_Box->Odd_D_On , &elem_num, Model->t, -1*sign, A_Basis);
   }
   
   //Need to fill zero in the diagnal elements for the inverse iteration method
   EXACT_MAKE_ELEM_ON(basis, 0, dim_onsite, Ham_Box->Zero_On, &elem_num, 0.0, A_Basis);

   for (i = 0; i < elem_num; i++) {
      A_Basis->Check[i] = -1;
   }
   
   A_Basis->elem_num = elem_num;
   
}
