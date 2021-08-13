//
//  Header.h
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/12.
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
//#include <mkl.h>

typedef struct {
   
   CRS1 **Ham         ;
   CRS1 **SpL_RE      ;
   CRS1 **SmL_RE      ;
   CRS1 **SzL_RE      ;
   CRS1 **CUp_1_RE    ;
   CRS1 **CDown_1_RE  ;
   CRS1 **CUp_1_D_RE  ;
   CRS1 **CDown_1_D_RE;
   CRS1 **CUp_2_RE    ;
   CRS1 **CDown_2_RE  ;
   CRS1 **CUp_2_D_RE  ;
   CRS1 **CDown_2_D_RE;
   
   CRS1 **SpL_LE    ;
   CRS1 **SmL_LE    ;
   CRS1 **SzL_LE    ;
   CRS1 **CUp_1_LE    ;
   CRS1 **CDown_1_LE  ;
   CRS1 **CUp_1_D_LE  ;
   CRS1 **CDown_1_D_LE;
   CRS1 **CUp_2_LE    ;
   CRS1 **CDown_2_LE  ;
   CRS1 **CUp_2_D_LE  ;
   CRS1 **CDown_2_D_LE;
   
   CRS1 *SzL_On       ;
   CRS1 *SpL_On       ;
   CRS1 *SmL_On       ;
   CRS1 *CUp_1_On     ;
   CRS1 *CDown_1_On   ;
   CRS1 *CUp_1_D_On   ;
   CRS1 *CDown_1_D_On ;
   CRS1 *CUp_2_On     ;
   CRS1 *CDown_2_On   ;
   CRS1 *CUp_2_D_On   ;
   CRS1 *CDown_2_D_On ;
   CRS1 *Ham_On       ;
   
   //For Expectation Values
   CRS1 *SzC_1_On;
   CRS1 *SzC_2_On;
   CRS1 *SxC_1_On;
   CRS1 *SxC_2_On;
   CRS1 *SxL_On  ;
   CRS1 *NC_1_On ;
   CRS1 *NC_2_On ;
   CRS1 *SC_1SL_On;
   CRS1 *SC_2SL_On;
   
   CCS1 **TM;
   CRS1 **TM_D;
   
   short **Basis_LL_LLLR;
   short **Basis_LR_LLLR;
   int ***Basis_Inv_LLLR;
   int *Dim_LLLR;
   
   int **Tot_Sz;
   int **Tot_Ele;
   int **Tot_Ele_1;
   int **Tot_Ele_2;
   int *Dim;
   
} BLOCK;

typedef struct {
   
   CRS1 *Ham_On      ;
   CRS1 *SpL_On      ;
   CRS1 *SmL_On      ;
   CRS1 *SzL_On      ;
   CRS1 *CUp_1_On    ;
   CRS1 *CUp_1_D_On  ;
   CRS1 *CDown_1_On  ;
   CRS1 *CDown_1_D_On;
   CRS1 *CUp_2_On    ;
   CRS1 *CUp_2_D_On  ;
   CRS1 *CDown_2_On  ;
   CRS1 *CDown_2_D_On;
   
   CRS1 *Ham_System         ;
   CRS1 *SpL_RE_System      ;
   CRS1 *SmL_RE_System      ;
   CRS1 *SzL_RE_System      ;
   CRS1 *CUp_1_RE_System    ;
   CRS1 *CDown_1_RE_System  ;
   CRS1 *CUp_1_D_RE_System  ;
   CRS1 *CDown_1_D_RE_System;
   CRS1 *CUp_2_RE_System    ;
   CRS1 *CDown_2_RE_System  ;
   CRS1 *CUp_2_D_RE_System  ;
   CRS1 *CDown_2_D_RE_System;
   
   CRS1 *SpL_LE_System      ;
   CRS1 *SmL_LE_System      ;
   CRS1 *SzL_LE_System      ;
   CRS1 *CUp_LE_System      ;
   CRS1 *CDown_1_LE_System  ;
   CRS1 *CUp_1_D_LE_System  ;
   CRS1 *CDown_1_D_LE_System;
   CRS1 *CUp_1_LE_System    ;
   CRS1 *CDown_2_LE_System  ;
   CRS1 *CUp_2_D_LE_System  ;
   CRS1 *CDown_2_D_LE_System;
   CRS1 *CUp_2_LE_System    ;
   
   CRS1 *Ham_Enviro       ;
   CRS1 *SpL_RE_Enviro    ;
   CRS1 *SmL_RE_Enviro    ;
   CRS1 *SzL_RE_Enviro    ;
   CRS1 *CUp_1_RE_Enviro    ;
   CRS1 *CDown_1_RE_Enviro  ;
   CRS1 *CUp_1_D_RE_Enviro  ;
   CRS1 *CDown_1_D_RE_Enviro;
   CRS1 *CUp_2_RE_Enviro    ;
   CRS1 *CDown_2_RE_Enviro  ;
   CRS1 *CUp_2_D_RE_Enviro  ;
   CRS1 *CDown_2_D_RE_Enviro;
   
   CRS1 *SpL_LE_Enviro    ;
   CRS1 *SmL_LE_Enviro    ;
   CRS1 *SzL_LE_Enviro    ;
   CRS1 *CUp_1_LE_Enviro    ;
   CRS1 *CDown_1_LE_Enviro  ;
   CRS1 *CUp_1_D_LE_Enviro  ;
   CRS1 *CDown_1_D_LE_Enviro;
   CRS1 *CUp_2_LE_Enviro    ;
   CRS1 *CDown_2_LE_Enviro  ;
   CRS1 *CUp_2_D_LE_Enviro  ;
   CRS1 *CDown_2_D_LE_Enviro;
   
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
   double h_z ;
   double mu  ;
   
   char *BC;
   
} HAM_BOX;


//EXPECTATIONS
void EXPECTATION_INTERSITE_Q2(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, int tot_sz, double **Temp_V1, double **Temp_V2, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status);
void EXPECTATION_SC_CORRELATIONS(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void EXPECTATION_VALUES(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status);
void STORE_T_MATRIX(DMRG_BASIS *Basis, DMRG_SYSTEM_INFO *Dmrg_System, BLOCK *System, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status, int LL_site, int dim_onsite);
void TRANSFORM_MATRIX_ONE_CORRELATIONS(CRS1 *M_On, CRS1 **Out1, CRS1 **Out2, BLOCK *System, MODEL_1DTKLM_VF *Model, char Enviro_Copy[100], char Name[100], int p_threads);
void TRANSFORM_MATRIX_TWO_CORRELATIONS(CRS1 *Onsite_Matrix, CRS1 **Out, BLOCK *System, MODEL_1DTKLM_VF *Model, int p_threads);
void SC_CORRELATIONS_ONSITE(int sc_sz, int sc_ele_1, int sc_ele_2, SC_MAT_1DTKLM_VF *Sc_Mat, CRS1 **M_Reference, CRS1 ***M_Enviro, MODEL_1DTKLM_VF *Model, DMRG_WHOLE_BASIS_Q4 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status, double *GS_Vec);

//STATUS
void OUTPUT_AVERAGE_VALUES(double *Val, char Name[], int start, int end, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_ENERGY(MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_FOURIER_COMPONENTS(double *Val, char Name[], int start, int end, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_INTERSITE_VALUES(double *Out, double *Onsite_Val, char Name[], MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_ONSITE_VALUES(double *Out, char Name[], MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_T_ERROR(DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DTKLM_VF *Model);
void OUTPUT_SC_CORRELATIONS(SC_MAT_1DTKLM_VF *Sc_Mat, int site_start, int site_end, char File_Name[], char D1[], char D2[], MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status);
void PRINT_STATUS(DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);

//ALLOCATIONS
void FREE_HAM_BOX(HAM_BOX **Box, MODEL_1DTKLM_VF *Model);
void FREE_T_MATRIX(BLOCK *System, BLOCK *Enviro, int tot_site, int sweep_now, int sweep, char Enviro_Copy[100]);
void FREE_BLOCK_MATRIX(BLOCK *Block, MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param);
void FREE_WHOLE_BASIS_Q2_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis, MODEL_1DTKLM_VF *Model);
void FREE_WHOLE_BASIS_Q4_SUPERBLOCK(DMRG_WHOLE_BASIS_Q4 *Basis, MODEL_1DTKLM_VF *Model);
void PRINT_MEM_STATUS(MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param);
DMRG_BASIS *GET_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
BLOCK *GET_BLOCK_MATRIX(MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param);
DMRG_STATUS *GET_DMRG_STATUS(MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param);
HAM_BOX **GET_HAM_BOX(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, int LL_site, int RR_site);
CRS1 *GET_HAM_LLLR(BLOCK *System, BLOCK *Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status, MODEL_1DTKLM_VF *Model);
CRS1 *GET_HAM_LLLRRRRL(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
DMRG_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
DMRG_WHOLE_BASIS_Q4 *GET_WHOLE_BASIS_Q4_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void GET_GROUND_STATE(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);

//MAIN
void DMRG(MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param);
int FIND_LLLRRRRL_DIM(int target_tot_ele_1, int target_tot_ele_2, int target_tot_sz, BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status);
double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DTKLM_VF *Model);
void MAKE_ELEMENT_HAM_LLLR(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box);
void MAKE_ELEMENT_HAM_LLLRRRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int *Inv, HAM_BOX *Box);
void MAKE_OP_EDGE(BLOCK *Block, MODEL_1DTKLM_VF *Model);
void RENORMALIZE(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void TRANSFORM_MATRIX(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DTKLM_VF *Model, DMRG_TIME *Dmrg_Time);




#endif /* Header_h */
