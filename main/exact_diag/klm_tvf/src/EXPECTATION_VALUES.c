//
//  EXPECTATION_VALUES.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/26.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_VALUES(MODEL_1DKLM_TVF *Model, EXACT_HAM_INFO *Ham_Info, EXACT_BASIS_INFO *Basis_Info, EXACT_TIME *Time) {
   
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

   EXACT_EXPECTATION_ONSITE(Ham_Box->SxL_On   , SxL   , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SxC_On   , SxC   , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SzLSzL_On, SzLSzL, Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SxLSxL_On, SxLSxL, Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->SzCSzC_On, SzCSzC, Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->NC_On    , NC    , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_ONSITE(Ham_Box->NCNC_On  , NCNC  , Ham_Info->Vector[0], Temp_V1, Model->dim_onsite, Model->tot_site, Model->p_threads, Basis_Info);
   
   OUTPUT_ONSITE_VALUES(SxL   , "SxL"   , Model);
   OUTPUT_ONSITE_VALUES(SxC   , "SxC"   , Model);
   OUTPUT_ONSITE_VALUES(SzLSzL, "SzLSzL", Model);
   OUTPUT_ONSITE_VALUES(SxLSxL, "SxLSxL", Model);
   OUTPUT_ONSITE_VALUES(SzCSzC, "SzCSzC", Model);
   OUTPUT_ONSITE_VALUES(NC    , "NC"    , Model);
   OUTPUT_ONSITE_VALUES(NCNC  , "NCNC"  , Model);
   
   OUTPUT_AVERAGE_VALUES(SxL   , "SxL"   , Model);
   OUTPUT_AVERAGE_VALUES(SxC   , "SxC"   , Model);
   OUTPUT_AVERAGE_VALUES(SzLSzL, "SzLSzL", Model);
   OUTPUT_AVERAGE_VALUES(SxLSxL, "SxLSxL", Model);
   OUTPUT_AVERAGE_VALUES(SzCSzC, "SzCSzC", Model);
   OUTPUT_AVERAGE_VALUES(NC    , "NC"    , Model);
   OUTPUT_AVERAGE_VALUES(NCNC  , "NCNC"  , Model);
   
   Time->exp_values = omp_get_wtime() - Time->exp_values;

   Time->cf = omp_get_wtime();
   
   EXACT_EXPECTATION_INTERSITE_Q0(Ham_Box->SxL_On, Ham_Box->SxL_On, Model->cf_origin, Model->tot_site, SxL_CF, Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_INTERSITE_Q0(Ham_Box->SxC_On, Ham_Box->SxC_On, Model->cf_origin, Model->tot_site, SxC_CF, Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, Basis_Info);
   EXACT_EXPECTATION_INTERSITE_Q0(Ham_Box->NC_On , Ham_Box->NC_On , Model->cf_origin, Model->tot_site, NC_CF , Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, Basis_Info);

   FREE_ARRAY_DOUBLE1(Temp_V1);
   FREE_ARRAY_DOUBLE1(Temp_V2);
   
   EXACT_WHOLE_BASIS_Q1 *W_Basis = GET_WHOLE_BASIS_Q1(Model, Time);
   Temp_V1 = GET_ARRAY_DOUBLE1(W_Basis->max_dim);
   Temp_V2 = GET_ARRAY_DOUBLE1(W_Basis->max_dim);
   EXPECTATION_INTERSITE_Q1(Ham_Box->SzL_On, Ham_Box->SzL_On, Model->cf_origin, Model->tot_site, SzL_CF, Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, W_Basis, Model->tot_parity);
   EXPECTATION_INTERSITE_Q1(Ham_Box->SzC_On, Ham_Box->SzC_On, Model->cf_origin, Model->tot_site, SzC_CF, Ham_Info->Vector[0], Temp_V1, Temp_V2, Model->dim_onsite, Model->p_threads, W_Basis, Model->tot_parity);
   
   OUTPUT_INTERSITE_VALUES(SxL_CF, SxL, Model->cf_origin, Model->tot_site, "SxL_CF", Model);
   OUTPUT_INTERSITE_VALUES(SxC_CF, SxC, Model->cf_origin, Model->tot_site, "SxC_CF", Model);
   OUTPUT_INTERSITE_VALUES(NC_CF , NC , Model->cf_origin, Model->tot_site, "NC_CF" , Model);
   OUTPUT_INTERSITE_VALUES(SzL_CF, SzL, Model->cf_origin, Model->tot_site, "SzL_CF", Model);
   OUTPUT_INTERSITE_VALUES(SzC_CF, SzC, Model->cf_origin, Model->tot_site, "SzC_CF", Model);
   
   FREE_WHOLE_BASIS_Q1(W_Basis);
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
   FREE_ARRAY_DOUBLE1(Temp_V1);
   FREE_ARRAY_DOUBLE1(Temp_V2);
   
   Time->cf = omp_get_wtime() - Time->cf;

}
