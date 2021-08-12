#include <stdio.h>
#include "SML.h"

void MATRIX_MATRIX_MATRIX_PRODUCT_CRS1(CRS1 *T_Matrix,
                                       CRS1 *Matrix,
                                       CRS1 *T_Dagger_Matrix,
                                       CRS1 *Out,
                                       CRS1 *Work
                                       ){
   
   MATRIX_PRODUCT_CRS1(Matrix, T_Matrix, Work);
   MATRIX_PRODUCT_CRS1(T_Dagger_Matrix, Work, Out);
   
}
