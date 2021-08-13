//
//  MAKE_ONSITE_HAM.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ONSITE_HAM(MODEL_1DXXZ_VF *Model, CRS1 *M) {
   
   int dim = Model->spin + 1;
   
   CRS1 *SzSz = GET_CRS1(dim, dim*dim);
   CRS1 *Sz   = GET_CRS1(dim, dim*dim);
   
   ONSITE_SZSZ_SZBASIS_HB(Model->spin, SzSz, 1.0);
   ONSITE_SZ_SZBASIS_HB(Model->spin, Sz, 1.0);

   int row;
   long iter;
   double D_z = Model->D_z;
   double h_z = Model->h_z;
   
   for (row = 0; row < dim; row++) {
      for (iter = SzSz->Row[row]; iter < SzSz->Row[row + 1]; iter++) {
         SzSz->Val[iter] = SzSz->Val[iter]*D_z;
      }
   }
   
   for (row = 0; row < dim; row++) {
      for (iter = Sz->Row[row]; iter < Sz->Row[row + 1]; iter++) {
         Sz->Val[iter] = Sz->Val[iter]*h_z;
      }
   }
   
   MATRIX_SUM_CRS1(SzSz, Sz, M);
   SORT_COLUMN_CRS1(M,1);
      
   FREE_CRS1(SzSz);
   FREE_CRS1(Sz);
}
