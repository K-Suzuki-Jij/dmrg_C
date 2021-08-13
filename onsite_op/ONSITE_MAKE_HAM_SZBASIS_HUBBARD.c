#include "model.h"
#include "SML.h"
#include "onsite.h"

void ONSITE_MAKE_HAM_SZBASIS_HUBBARD(MODEL_1DHUBBARD_VF *Model, CRS1 *M) {
   
   int dim      = Model->dim_onsite;
   
   CRS1 *SzC_hz    = GET_CRS1(dim, dim*dim);
   CRS1 *NCNC      = GET_CRS1(dim, dim*dim);
   CRS1 *N_mu      = GET_CRS1(dim, dim*dim);
   CRS1 *M1        = GET_CRS1(dim, dim*dim);

   ONSITE_SZC_SZBASIS_HUBBARD          (SzC_hz, Model->h_z);
   ONSITE_NC_UP_NC_DOWN_SZBASIS_HUBBARD(NCNC  , Model->U  );
   ONSITE_NC_SZBASIS_HUBBARD           (N_mu  , Model->mu );
   
   MATRIX_SUM_CRS1(SzC_hz  , NCNC, M1);
   MATRIX_SUM_CRS1(M1, N_mu  , M);

   FREE_CRS1(SzC_hz);
   FREE_CRS1(NCNC  );
   FREE_CRS1(N_mu  );
   FREE_CRS1(M1    );

}
