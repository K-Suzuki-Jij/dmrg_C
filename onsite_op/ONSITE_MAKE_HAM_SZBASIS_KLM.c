#include "model.h"
#include "SML.h"
#include "onsite.h"

void ONSITE_MAKE_HAM_SZBASIS_KLM(MODEL_1DKLM_VF *Model, CRS1 *M) {
   
   int dim      = Model->dim_onsite;
   int spin_loc = Model->spin_loc;
   
   CRS1 *SCSL_J    = GET_CRS1(dim, dim*dim);
   CRS1 *SzC_hz    = GET_CRS1(dim, dim*dim);
   CRS1 *SzL_hz    = GET_CRS1(dim, dim*dim);
   CRS1 *SzLSzL_Dz = GET_CRS1(dim, dim*dim);
   CRS1 *N_mu      = GET_CRS1(dim, dim*dim);
   CRS1 *M1        = GET_CRS1(dim, dim*dim);
   CRS1 *M2        = GET_CRS1(dim, dim*dim);
   CRS1 *M3        = GET_CRS1(dim, dim*dim);

   ONSITE_SCSL_SZBASIS_KLM  (spin_loc, SCSL_J   , Model->J);
   ONSITE_SZC_SZBASIS_KLM   (spin_loc, SzC_hz   , Model->h_z);
   ONSITE_SZL_SZBASIS_KLM   (spin_loc, SzL_hz   , Model->h_z);
   ONSITE_SZLSZL_SZBASIS_KLM(spin_loc, SzLSzL_Dz, Model->D_z);
   ONSITE_NC_SZBASIS_KLM    (spin_loc, N_mu     , Model->mu);

   
   MATRIX_SUM_CRS1(SCSL_J  , SzC_hz   , M1);
   MATRIX_SUM_CRS1(SzL_hz  , SzLSzL_Dz, M2);
   MATRIX_SUM_CRS1(M1, M2  , M3);
   MATRIX_SUM_CRS1(M3, N_mu, M);

   FREE_CRS1(SCSL_J   );
   FREE_CRS1(SzC_hz   );
   FREE_CRS1(SzL_hz   );
   FREE_CRS1(SzLSzL_Dz);
   FREE_CRS1(N_mu     );
   FREE_CRS1(M1       );
   FREE_CRS1(M2       );
   FREE_CRS1(M3       );

}
