//
//  Header.h
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/13.
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
#include <assert.h>
#include "dmrg.h"
#include "SML.h"
#include "onsite.h"
#include "model.h"

typedef struct {
   
   CRS1 **Ham      ;
   CRS1 **SpL_RE   ;
   CRS1 **SmL_RE   ;
   CRS1 **SzL_RE   ;
   CRS1 **Even_RE  ;
   CRS1 **Odd_RE   ;
   CRS1 **Even_D_RE;
   CRS1 **Odd_D_RE ;
   
   CRS1 **SpL_LE   ;
   CRS1 **SmL_LE   ;
   CRS1 **SzL_LE   ;
   CRS1 **Even_LE  ;
   CRS1 **Odd_LE   ;
   CRS1 **Even_D_LE;
   CRS1 **Odd_D_LE ;
   
   CRS1 *SzL_On   ;
   CRS1 *SpL_On   ;
   CRS1 *SmL_On   ;
   CRS1 *Even_On  ;
   CRS1 *Odd_On   ;
   CRS1 *Even_D_On;
   CRS1 *Odd_D_On ;
   CRS1 *Ham_On   ;
   
   //For Expectation Values
   CRS1 *SxL_On;
   CRS1 *SzC_On;
   CRS1 *SxC_On;
   CRS1 *SCSL_On;
   CRS1 *NC_On;
   
   CCS1 **TM;
   CRS1 **TM_D;
   
   short **Basis_LL_LLLR;
   short **Basis_LR_LLLR;
   int ***Basis_Inv_LLLR;
   int *Dim_LLLR;
   
   int **Tot_Parity;
   int **Tot_Ele;
   int *Dim;
   
} BLOCK;

typedef struct {
   
   CRS1 *Ham_On   ;
   CRS1 *SpL_On   ;
   CRS1 *SmL_On   ;
   CRS1 *SzL_On   ;
   CRS1 *Even_On  ;
   CRS1 *Even_D_On;
   CRS1 *Odd_On   ;
   CRS1 *Odd_D_On ;
   
   CRS1 *Ham_System      ;
   CRS1 *SpL_RE_System   ;
   CRS1 *SmL_RE_System   ;
   CRS1 *SzL_RE_System   ;
   CRS1 *Even_RE_System  ;
   CRS1 *Odd_RE_System   ;
   CRS1 *Even_D_RE_System;
   CRS1 *Odd_D_RE_System ;
   
   CRS1 *SpL_LE_System   ;
   CRS1 *SmL_LE_System   ;
   CRS1 *SzL_LE_System   ;
   CRS1 *Even_LE_System  ;
   CRS1 *Odd_LE_System   ;
   CRS1 *Even_D_LE_System;
   CRS1 *Odd_D_LE_System ;
   
   CRS1 *Ham_Enviro      ;
   CRS1 *SpL_RE_Enviro   ;
   CRS1 *SmL_RE_Enviro   ;
   CRS1 *SzL_RE_Enviro   ;
   CRS1 *Even_RE_Enviro  ;
   CRS1 *Odd_RE_Enviro   ;
   CRS1 *Even_D_RE_Enviro;
   CRS1 *Odd_D_RE_Enviro ;
   
   CRS1 *SpL_LE_Enviro   ;
   CRS1 *SmL_LE_Enviro   ;
   CRS1 *SzL_LE_Enviro   ;
   CRS1 *Even_LE_Enviro  ;
   CRS1 *Odd_LE_Enviro   ;
   CRS1 *Even_D_LE_Enviro;
   CRS1 *Odd_D_LE_Enviro ;
   
   int *Ele_LL;
   int *Ele_RR;
   int *Ele_On;
   
   double SSD_LR;
   double SSD_RL;
   double SSD_LLLR;
   double SSD_LRRL;
   double SSD_RRRL;
   
   double t   ;
   double J   ;
   double D_z ;
   double I_xy;
   double I_z ;
   double h_xc;
   double h_xl;
   double mu  ;
   
   char *BC;   
   
} HAM_BOX;

//STATUS
void OUTPUT_ENERGY(MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_ONSITE_VALUES(double *Out, char Name[], MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_AVERAGE_VALUES(double *Val, char Name[], int start, int end, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_FOURIER_COMPONENTS(double *Val, char Name[], int start, int end, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_INTERSITE_VALUES(double *Out, double *Onsite_Val, char Name[], MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_SC_CORRELATIONS(SC_MAT_1DKLM_TVF *Sc_Mat, int site_start, int site_end, char File_Name[], char D1[], char D2[], MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_T_ERROR(DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DKLM_TVF *Model);
void PRINT_STATUS(DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void PRINT_MEM_STATUS(MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param);

//ALLOCATIONS
void FREE_BLOCK_MATRIX(BLOCK *Block, MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param);
void FREE_HAM_BOX(HAM_BOX **Box, MODEL_1DKLM_TVF *Model);
void FREE_WHOLE_BASIS_Q2_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis);
void FREE_WHOLE_BASIS_Q1_SUPERBLOCK(DMRG_WHOLE_BASIS_Q1 *Basis);
CRS1 *GET_HAM_LLLR(BLOCK *Block_System, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DKLM_TVF *Model);
DMRG_BASIS *GET_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
CRS1 *GET_HAM_LLLRRRRL(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
HAM_BOX **GET_HAM_BOX(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, int LL_site, int RR_site);
DMRG_STATUS *GET_DMRG_STATUS(MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param);
BLOCK *GET_BLOCK_MATRIX(MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param);
void GET_GROUND_STATE(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void FREE_T_MATRIX(BLOCK *System, BLOCK *Enviro, int tot_site, int sweep_now, int sweep, char Enviro_Copy[100]);

//MAIN
void DMRG(MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param);
int FIND_LLLRRRRL_DIM(int target_tot_ele, int target_tot_parity, BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status);
double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DKLM_TVF *Model);
void MAKE_ELEMENT_HAM_LLLR(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box);
void MAKE_ELEMENT_HAM_LLLRRRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int *Inv, HAM_BOX *Box);
void MAKE_OP_EDGE(BLOCK *Block, MODEL_1DKLM_TVF *Model);
void RENORMALIZE(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void TRANSFORM_MATRIX(BLOCK *Block_System, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DKLM_TVF *Model, DMRG_TIME *Dmrg_Time);
DMRG_WHOLE_BASIS_Q1 *GET_WHOLE_BASIS_Q1_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
DMRG_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);

//EXPECTATIONS
void EXPECTATION_INTERSITE_PARITY(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, int tot_parity, double *Temp_V1, double *Temp_V2, int p_threads, DMRG_WHOLE_BASIS_Q1 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status);
void EXPECTATION_SC_CORRELATIONS(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void EXPECTATION_VALUES(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status);
void SC_CORRELATIONS(int sc_p, char LS_Couple[100], BLOCK *System, BLOCK *Enviro, SC_MAT_1DKLM_TVF *Sc_Mat, MODEL_1DKLM_TVF *Model, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status, double *GS_Vec);
void TRANSFORM_MATRIX_FOR_EXPECTATION_VALUES(BLOCK *Block_System, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void TRANSFORM_MATRIX_FOR_SC_CORRELATIONS(BLOCK *Block_System, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void STORE_T_MATRIX(DMRG_BASIS *Basis, DMRG_SYSTEM_INFO *Dmrg_System, BLOCK *System, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status, int LL_site, int dim_onsite);
void GREENS_FUNCTION(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status);


#endif /* Header_h */
