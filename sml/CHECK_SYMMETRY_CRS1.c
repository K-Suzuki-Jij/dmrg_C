//
//  CHECK_SYMMETRY_CRS1_2.c
//  Diag
//
//  Created by Kohei Suzuki on 2018/12/04.
//  Copyright © 2018 Kohei Suzuki. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "SML.h"

int CHECK_SYMMETRY_CRS1(CRS1 *M, double zero, int p_threads) {
   
   long i,j,inv;
   
#pragma omp parallel for private (j,inv) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      for (j = M->Row[i]; j < M->Row[i+1]; j++) {
         inv = BINARY_SEARCH_INT1(M->Col, M->Row[M->Col[j]], M->Row[M->Col[j] + 1], (int)i);
         if (inv == -1) {
            printf("The input matrix is not symmetric.\n");
            printf("Corresponding element does not exist\n");
            printf("row=%ld,col=%d,val=%.15lf\n",i,M->Col[j],M->Val[j]);
            exit(1);
         }
         if (fabs(M->Val[j] - M->Val[inv]) > zero) {
            printf("The input matrix is not symmetric.\n");
            printf("M[%d][%d]=%.15lf,%.15lf=M[%d][%d]\n",
                   (int)i,M->Col[j],M->Val[j],
                   M->Val[inv],M->Col[j],(int)i
                   );
            exit(1);
         }
      }
   }
   
   return 0;
   
}
