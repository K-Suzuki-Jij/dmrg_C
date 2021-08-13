#include "model.h"
#include "SML.h"
#include "onsite.h"

void ONSITE_MAKE_HAM_SPIN_INV_BASIS_KLM(MODEL_1DKLM_TVF *Model, CRS1 *M) {
   
   int dim      = Model->dim_onsite;
   int spin_loc = Model->spin_loc;
   
   CRS1 *SCSL_J    = GET_CRS1(dim, dim*dim);
   CRS1 *SxC_hx    = GET_CRS1(dim, dim*dim);
   CRS1 *SxL_hx    = GET_CRS1(dim, dim*dim);
   CRS1 *SzLSzL_Dz = GET_CRS1(dim, dim*dim);
   CRS1 *N_mu      = GET_CRS1(dim, dim*dim);
   CRS1 *M1        = GET_CRS1(dim, dim*dim);
   CRS1 *M2        = GET_CRS1(dim, dim*dim);
   CRS1 *M3        = GET_CRS1(dim, dim*dim);

   ONSITE_SCSL_SPIN_INV_BASIS_KLM  (spin_loc, SCSL_J   , Model->J);
   ONSITE_SXC_SPIN_INV_BASIS_KLM   (spin_loc, SxC_hx   , Model->h_xc);
   ONSITE_SXL_SPIN_INV_BASIS_KLM   (spin_loc, SxL_hx   , Model->h_xl);
   ONSITE_SZLSZL_SPIN_INV_BASIS_KLM(spin_loc, SzLSzL_Dz, Model->D_z);
   ONSITE_NC_SPIN_INV_BASIS_KLM    (spin_loc, N_mu     , Model->mu);

   MATRIX_SUM_CRS1(SCSL_J  , SxC_hx   , M1);
   MATRIX_SUM_CRS1(SxL_hx  , SzLSzL_Dz, M2);
   MATRIX_SUM_CRS1(M1, M2  , M3);
   MATRIX_SUM_CRS1(M3, N_mu, M);
   
   FREE_CRS1(SCSL_J   );
   FREE_CRS1(SxC_hx   );
   FREE_CRS1(SxL_hx   );
   FREE_CRS1(SzLSzL_Dz);
   FREE_CRS1(N_mu     );
   FREE_CRS1(M1       );
   FREE_CRS1(M2       );
   FREE_CRS1(M3       );

}
