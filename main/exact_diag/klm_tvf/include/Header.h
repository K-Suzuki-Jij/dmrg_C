//
//  Header.h
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/25.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#ifndef Header_h
#define Header_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <limits.h>
//#include <mkl.h>
#include <assert.h>
#include "exact.h"
#include "SML.h"
#include "onsite.h"
#include "model.h"

typedef struct {
   
   CRS1 *Ham_On   ;
   CRS1 *SpL_On   ;
   CRS1 *SmL_On   ;
   CRS1 *SzL_On   ;
   CRS1 *Even_On  ;
   CRS1 *Even_D_On;
   CRS1 *Odd_On   ;
   CRS1 *Odd_D_On ;
   CRS1 *Zero_On  ;
   
   CRS1 *SxL_On;
   CRS1 *SxC_On;
   CRS1 *SzC_On;
   CRS1 *NC_On;
   CRS1 *NCNC_On;
   CRS1 *SzLSzL_On;
   CRS1 *SxLSxL_On;
   CRS1 *SzCSzC_On;
   
   CRS1 **CCSL_On;
   CRS1 **CSL_On;
   
} HAM_BOX;


int FIND_DIM(MODEL_1DKLM_TVF *Model, EXACT_TIME *Time);
long *GET_BASIS(MODEL_1DKLM_TVF *Model, EXACT_TIME *Time, int dim);
void EXACT_DIAGONALIZATION(MODEL_1DKLM_TVF *Model, EXACT_PARAMETER *Param, EXACT_BASIS_INFO *Basis_Info, EXACT_HAM_INFO *Ham_Info, EXACT_TIME *Time);
CRS1 *GET_HAM(MODEL_1DKLM_TVF *Model, EXACT_BASIS_INFO *Basis_Info, EXACT_PARAMETER *Param, EXACT_TIME *Time);
HAM_BOX *GET_HAM_BOX(MODEL_1DKLM_TVF *Model);
void MAKE_ELEMENT_HAM(long basis, EXACT_A_BASIS *A_Basis, HAM_BOX *Ham_Box, MODEL_1DKLM_TVF *Model);
void FREE_HAM_BOX(HAM_BOX *Box, MODEL_1DKLM_TVF *Model);
void EXPECTATION_VALUES(MODEL_1DKLM_TVF *Model, EXACT_HAM_INFO *Ham_Info, EXACT_BASIS_INFO *Basis_Info, EXACT_TIME *Time);
EXACT_WHOLE_BASIS_Q1 *GET_WHOLE_BASIS_Q1(MODEL_1DKLM_TVF *Model, EXACT_TIME *Time);
void EXPECTATION_INTERSITE_Q1(CRS1 *M_O, CRS1 *M_R, int start, int end, double *Out, double *Vec, double *T_Vec1, double *T_Vec2, int dim_onsite, int p_threads, EXACT_WHOLE_BASIS_Q1 *W_Basis, int tot_parity);
void FREE_WHOLE_BASIS_Q1(EXACT_WHOLE_BASIS_Q1 *W_Basis);
void OUTPUT_ONSITE_VALUES(double *Out, char Name[], MODEL_1DKLM_TVF *Model);
void OUTPUT_INTERSITE_VALUES(double *Out, double *Onsite_Val, int origin, int end, char Name[], MODEL_1DKLM_TVF *Model);
void OUTPUT_ENERGY(MODEL_1DKLM_TVF *Model, EXACT_HAM_INFO *Ham_Info);
void PRINT_STATUS(MODEL_1DKLM_TVF *Model, EXACT_TIME *Time, EXACT_HAM_INFO *Ham_Info, EXACT_PARAMETER *Param);
void OUTPUT_AVERAGE_VALUES(double *Out, char Name[], MODEL_1DKLM_TVF *Model);
EXACT_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2(MODEL_1DKLM_TVF *Model, EXACT_TIME *Time);
void EXPECTATION_SC_CORRELATIONS(MODEL_1DKLM_TVF *Model, EXACT_HAM_INFO *Ham_Info, EXACT_BASIS_INFO *Basis_Info, EXACT_TIME *Time);
void SC_CORRELATIONS(int sc_p, int site_start, double *GS_Vec, SC_MAT_1DKLM_TVF *Sc_Mat, EXACT_WHOLE_BASIS_Q2 *W_Basis, MODEL_1DKLM_TVF *Model);
void OUTPUT_SC_CORRELATIONS(SC_MAT_1DKLM_TVF *Sc_Mat, int site_start, int site_end, char File_Name[], char D1[], char D2[], MODEL_1DKLM_TVF *Model);
void FREE_WHOLE_BASIS_Q2(EXACT_WHOLE_BASIS_Q2 *W_Basis);

#endif /* Header_h */
