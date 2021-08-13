//
//  EXPECTATION_VALUES.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/11.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_VALUES(MODEL_1DKLM_VF *Model, EXACT_HAM_INFO *Ham_Info, EXACT_BASIS_INFO *Basis_Info, EXACT_TIME *Time) {
   
   Time->exp_values = omp_get_wtime();
   
   HAM_BOX *Ham_Box = GET_HAM_BOX(Model);
   double *SxC      = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxL      = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzC      = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzL      = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzLSzL   = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxLSxL   = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzCSzC   = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *NC       = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *NCNC     = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzL_CF   = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzC_CF   = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxL_CF   = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxC_CF   = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *NC_CF    = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *Temp_V1  = GET_ARRAY_DOUBLE1(Basis_Info->dim);
   double *Temp_V2  = GET_ARRAY_DOUBLE1(Basis_Info->dim);
   
   EXACT_EXPECTATION_ONSITE(Ham_Box->SzL_On   , SzL   , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SzC_On   , SzC   , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SzLSzL_On, SzLSzL, Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SxLSxL_On, SxLSxL, Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SzCSzC_On, SzCSzC, Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->NC_On    , NC    , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->NCNC_On  , NCNC  , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);

   OUTPUT_ONSITE_VALUES(SzL   , "SzL"   , Model);
   OUTPUT_ONSITE_VALUES(SzC   , "SzC"   , Model);
   OUTPUT_ONSITE_VALUES(SzLSzL, "SzLSzL", Model);
   OUTPUT_ONSITE_VALUES(SxLSxL, "SxLSxL", Model);
   OUTPUT_ONSITE_VALUES(SzCSzC, "SzCSzC", Model);
   OUTPUT_ONSITE_VALUES(NC    , "NC"    , Model);
   OUTPUT_ONSITE_VALUES(NCNC  , "NCNC"  , Model);
   
   OUTPUT_AVERAGE_VALUES(SzL   , "SzL"   , Model);
   OUTPUT_AVERAGE_VALUES(SzC   , "SzC"   , Model);
   OUTPUT_AVERAGE_VALUES(SzLSzL, "SzLSzL", Model);
   OUTPUT_AVERAGE_VALUES(SxLSxL, "SxLSxL", Model);
   OUTPUT_AVERAGE_VALUES(SzCSzC, "SzCSzC", Model);
   OUTPUT_AVERAGE_VALUES(NC    , "NC"    , Model);
   OUTPUT_AVERAGE_VALUES(NCNC  , "NCNC"  , Model);
   
   Time->exp_values = omp_get_wtime() - Time->exp_values;
   
   Time->cf = omp_get_wtime();
   
   EXACT_EXPECTATION_INTERSITE_Q0(Ham_Box->SzL_On, Ham_Box->SzL_On, Model->cf_origin, Model->tot_site, SzL_CF, Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_INTERSITE_Q0(Ham_Box->SzC_On, Ham_Box->SzC_On, Model->cf_origin, Model->tot_site, SzC_CF, Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_INTERSITE_Q0(Ham_Box->NC_On , Ham_Box->NC_On , Model->cf_origin, Model->tot_site, NC_CF , Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, Basis_Info);

   OUTPUT_INTERSITE_VALUES(SzC_CF, SzC, Model->cf_origin, Model->tot_site, "SzC_CF", Model);
   OUTPUT_INTERSITE_VALUES(SzL_CF, SzL, Model->cf_origin, Model->tot_site, "SzL_CF", Model);
   OUTPUT_INTERSITE_VALUES(NC_CF , NC , Model->cf_origin, Model->tot_site, "NC_CF" , Model);
   
   FREE_ARRAY_DOUBLE1(Temp_V1);
   FREE_ARRAY_DOUBLE1(Temp_V2);
   
   EXACT_WHOLE_BASIS_Q2 *W_Basis = GET_WHOLE_BASIS_Q2(Model, Time);
   double **Temp_VV1  = GET_ARRAY_DOUBLE2(2, W_Basis->max_dim);
   double **Temp_VV2  = GET_ARRAY_DOUBLE2(2, W_Basis->max_dim);
   EXPECTATION_INTERSITE_Q2(Ham_Box->SxL_On, Ham_Box->SxL_On, Model->cf_origin, Model->tot_site, SxL_CF, Ham_Info->Vector[0], Temp_VV1, Temp_VV2, Model->dim_onsite, Model->p_threads, W_Basis, Model->tot_sz);
   EXPECTATION_INTERSITE_Q2(Ham_Box->SxC_On, Ham_Box->SxC_On, Model->cf_origin, Model->tot_site, SxC_CF, Ham_Info->Vector[0], Temp_VV1, Temp_VV2, Model->dim_onsite, Model->p_threads, W_Basis, Model->tot_sz);
   
   OUTPUT_INTERSITE_VALUES(SxC_CF, SxC, Model->cf_origin, Model->tot_site, "SxC_CF", Model);
   OUTPUT_INTERSITE_VALUES(SxL_CF, SxL, Model->cf_origin, Model->tot_site, "SxL_CF", Model);
   
   FREE_WHOLE_BASIS_Q2(W_Basis, Model->tot_site*Model->spin_loc + Model->tot_ele);
   FREE_HAM_BOX(Ham_Box, Model);
   FREE_ARRAY_DOUBLE1(SxC    );
   FREE_ARRAY_DOUBLE1(SxL    );
   FREE_ARRAY_DOUBLE1(SzC    );
   FREE_ARRAY_DOUBLE1(SzL    );
   FREE_ARRAY_DOUBLE1(SzLSzL );
   FREE_ARRAY_DOUBLE1(SxLSxL );
   FREE_ARRAY_DOUBLE1(SzCSzC );
   FREE_ARRAY_DOUBLE1(NC     );
   FREE_ARRAY_DOUBLE1(NCNC   );
   FREE_ARRAY_DOUBLE1(SzL_CF );
   FREE_ARRAY_DOUBLE1(SzC_CF );
   FREE_ARRAY_DOUBLE1(SxL_CF );
   FREE_ARRAY_DOUBLE1(SxC_CF );
   FREE_ARRAY_DOUBLE1(NC_CF  );
   FREE_ARRAY_DOUBLE2(Temp_VV1, 2);
   FREE_ARRAY_DOUBLE2(Temp_VV2, 2);
   
   Time->cf = omp_get_wtime() - Time->cf;
   
}
