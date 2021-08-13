#include "model.h"
#include "SML.h"
#include "onsite.h"

void ONSITE_MAKE_HAM_SZBASIS_TKLM(MODEL_1DTKLM_VF *Model, CRS1 *M) {
   
   int dim      = Model->dim_onsite;
   int spin_loc = Model->spin_loc;
   
   CRS1 *SC_1SL_J  = GET_CRS1(dim, dim*dim);
   CRS1 *SC_2SL_J  = GET_CRS1(dim, dim*dim);
   CRS1 *SzC_1_hz  = GET_CRS1(dim, dim*dim);
   CRS1 *SzC_2_hz  = GET_CRS1(dim, dim*dim);
   CRS1 *SzL_hz    = GET_CRS1(dim, dim*dim);
   CRS1 *SzLSzL_Dz = GET_CRS1(dim, dim*dim);
   CRS1 *NC_1_mu   = GET_CRS1(dim, dim*dim);
   CRS1 *NC_2_mu   = GET_CRS1(dim, dim*dim);

   CRS1 *M1        = GET_CRS1(dim, dim*dim);
   CRS1 *M2        = GET_CRS1(dim, dim*dim);
   CRS1 *M3        = GET_CRS1(dim, dim*dim);
   CRS1 *M4        = GET_CRS1(dim, dim*dim);
   CRS1 *MM1       = GET_CRS1(dim, dim*dim);
   CRS1 *MM2       = GET_CRS1(dim, dim*dim);

   ONSITE_SC_1SL_SZBASIS_TKLM(spin_loc, SC_1SL_J  , Model->J);
   ONSITE_SC_2SL_SZBASIS_TKLM(spin_loc, SC_2SL_J  , Model->J);
   ONSITE_SZC_1_SZBASIS_TKLM (spin_loc, SzC_1_hz  , Model->h_z);
   ONSITE_SZC_2_SZBASIS_TKLM (spin_loc, SzC_2_hz  , Model->h_z);
   ONSITE_SZL_SZBASIS_TKLM   (spin_loc, SzL_hz    , Model->h_z);
   ONSITE_SZLSZL_SZBASIS_TKLM(spin_loc, SzLSzL_Dz , Model->D_z);
   ONSITE_NC_1_SZBASIS_TKLM  (spin_loc, NC_1_mu   , Model->mu);
   ONSITE_NC_2_SZBASIS_TKLM  (spin_loc, NC_2_mu   , Model->mu);

   MATRIX_SUM_CRS1(SC_1SL_J, SC_2SL_J  , M1);
   MATRIX_SUM_CRS1(SzC_1_hz, SzC_2_hz  , M2);
   MATRIX_SUM_CRS1(SzL_hz  , SzLSzL_Dz , M3);
   MATRIX_SUM_CRS1(NC_1_mu , NC_1_mu   , M4);
   
   MATRIX_SUM_CRS1(M1 , M2  , MM1);
   MATRIX_SUM_CRS1(M3 , M4  , MM2);
   MATRIX_SUM_CRS1(MM1, MM2 , M);
   
   FREE_CRS1(SC_1SL_J );
   FREE_CRS1(SC_2SL_J );
   FREE_CRS1(SzC_1_hz );
   FREE_CRS1(SzC_2_hz );
   FREE_CRS1(SzL_hz   );
   FREE_CRS1(SzLSzL_Dz);
   FREE_CRS1(NC_1_mu  );
   FREE_CRS1(NC_2_mu  );
   FREE_CRS1(M1       );
   FREE_CRS1(M2       );
   FREE_CRS1(M3       );
   FREE_CRS1(M4       );
   FREE_CRS1(MM1      );
   FREE_CRS1(MM2      );


}
