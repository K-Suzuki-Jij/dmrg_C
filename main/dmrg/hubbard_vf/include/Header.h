//
//  Header.h
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
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
   
   CRS1 **Ham         ;
   CRS1 **CUp_RE      ;
   CRS1 **CDown_RE    ;
   CRS1 **CUp_D_RE    ;
   CRS1 **CDown_D_RE  ;
   CRS1 **NC_Up_RE    ;
   CRS1 **NC_Down_RE  ;

   CRS1 **CUp_LE      ;
   CRS1 **CDown_LE    ;
   CRS1 **CUp_D_LE    ;
   CRS1 **CDown_D_LE  ;
   CRS1 **NC_Up_LE    ;
   CRS1 **NC_Down_LE  ;
   
   CRS1 *CUp_On     ;
   CRS1 *CDown_On   ;
   CRS1 *CUp_D_On   ;
   CRS1 *CDown_D_On ;
   CRS1 *NC_Up_On   ;
   CRS1 *NC_Down_On ;
   CRS1 *Ham_On     ;
   
   //For Expectation Values
   CRS1 *SzC_On;
   CRS1 *SxC_On;
   CRS1 *NC_On;
   CRS1 *DO_On;

   CCS1 **TM;
   CRS1 **TM_D;
   
   short **Basis_LL_LLLR;
   short **Basis_LR_LLLR;
   int ***Basis_Inv_LLLR;
   int *Dim_LLLR;
   
   int **Tot_Sz ;
   int **Tot_Ele;
   int *Dim;
   
} BLOCK;

typedef struct {
   
   CRS1 *Ham_On    ;
   CRS1 *CUp_On    ;
   CRS1 *CUp_D_On  ;
   CRS1 *CDown_On  ;
   CRS1 *CDown_D_On;
   CRS1 *NC_Up_On  ;
   CRS1 *NC_Down_On;

   CRS1 *Ham_System         ;
   CRS1 *CUp_RE_System      ;
   CRS1 *CDown_RE_System    ;
   CRS1 *CUp_D_RE_System    ;
   CRS1 *CDown_D_RE_System  ;
   CRS1 *NC_Up_RE_System    ;
   CRS1 *NC_Down_RE_System  ;
   
   CRS1 *CUp_LE_System      ;
   CRS1 *CDown_LE_System    ;
   CRS1 *CUp_D_LE_System    ;
   CRS1 *CDown_D_LE_System  ;
   CRS1 *NC_Up_LE_System    ;
   CRS1 *NC_Down_LE_System  ;
   
   CRS1 *Ham_Enviro         ;
   CRS1 *CUp_RE_Enviro      ;
   CRS1 *CDown_RE_Enviro    ;
   CRS1 *CUp_D_RE_Enviro    ;
   CRS1 *CDown_D_RE_Enviro  ;
   CRS1 *NC_Up_RE_Enviro    ;
   CRS1 *NC_Down_RE_Enviro  ;
   
   CRS1 *SpL_LE_Enviro      ;
   CRS1 *SmL_LE_Enviro      ;
   CRS1 *SzL_LE_Enviro      ;
   CRS1 *CUp_LE_Enviro      ;
   CRS1 *CDown_LE_Enviro    ;
   CRS1 *CUp_D_LE_Enviro    ;
   CRS1 *CDown_D_LE_Enviro  ;
   CRS1 *NC_Up_LE_Enviro    ;
   CRS1 *NC_Down_LE_Enviro  ;
   
   int *Ele_LL;
   int *Ele_RR;
   int *Ele_On;
   
   double SSD_LR;
   double SSD_RL;
   double SSD_LLLR;
   double SSD_LRRL;
   double SSD_RRRL;
   double SSD_LRRR;

   double t1 ;
   double t2 ;
   double U  ;
   double V  ;
   double h_z;
   double mu ;
   
   char *BC;
   
} HAM_BOX;



void DMRG(MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param);
BLOCK *GET_BLOCK_MATRIX(MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param);
DMRG_STATUS *GET_DMRG_STATUS(MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param);
void MAKE_OP_EDGE(BLOCK *Block, MODEL_1DHUBBARD_VF *Model);
DMRG_BASIS *GET_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
int FIND_LLLRRRRL_DIM(int target_tot_ele, int target_tot_sz, BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status);
void RENORMALIZE(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void GET_GROUND_STATE(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
HAM_BOX **GET_HAM_BOX(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, int LL_site, int RR_site);
CRS1 *GET_HAM_LLLRRRRL(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void MAKE_ELEMENT_HAM_LLLRRRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int *Inv, HAM_BOX *Box);
void FREE_HAM_BOX(HAM_BOX **Box, MODEL_1DHUBBARD_VF *Model);
void OUTPUT_T_ERROR(DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DHUBBARD_VF *Model);
void TRANSFORM_MATRIX(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DHUBBARD_VF *Model, DMRG_TIME *Dmrg_Time);
CRS1 *GET_HAM_LLLR(BLOCK *System, BLOCK *Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status, MODEL_1DHUBBARD_VF *Model);
double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DHUBBARD_VF *Model);
void MAKE_ELEMENT_HAM_LLLR(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box);
void STORE_T_MATRIX(DMRG_BASIS *Basis, DMRG_SYSTEM_INFO *Dmrg_System, BLOCK *System, MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status, int LL_site, int dim_onsite);
void PRINT_STATUS(DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void OUTPUT_ENERGY(MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_ONSITE_VALUES(double *Out, char Name[], MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_FOURIER_COMPONENTS(double *Val, char Name[], int start, int end, MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_AVERAGE_VALUES(double *Val, char Name[], int start, int end, MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status);
DMRG_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void FREE_WHOLE_BASIS_Q2_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis, MODEL_1DHUBBARD_VF *Model);
void EXPECTATION_INTERSITE_Q2(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, int tot_sz, double **Temp_V1, double **Temp_V2, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status);
void OUTPUT_INTERSITE_VALUES(double *Out, double *Onsite_Val, char Name[], MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status);
void EXPECTATION_VALUES(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status);
void FREE_BLOCK_MATRIX(BLOCK *Block, MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param);
void EXPECTATION_SC_CORRELATIONS(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
DMRG_WHOLE_BASIS_Q3 *GET_WHOLE_BASIS_Q3_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void FREE_WHOLE_BASIS_Q3_SUPERBLOCK(DMRG_WHOLE_BASIS_Q3 *Basis, MODEL_1DHUBBARD_VF *Model);
void SC_CORRELATIONS(int sc_sz, BLOCK *System, BLOCK *Enviro, SC_MAT_1DHUBBARD_VF *Sc_Mat, MODEL_1DHUBBARD_VF *Model, DMRG_WHOLE_BASIS_Q3 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status, double *GS_Vec);
void OUTPUT_SC_CORRELATIONS(SC_MAT_1DHUBBARD_VF *Sc_Mat, int site_start, int site_end, char File_Name[], char D1[], char D2[], MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status);
void FREE_T_MATRIX(BLOCK *System, BLOCK *Enviro, int tot_site, int sweep_now, int sweep, char Enviro_Copy[100]);
#endif /* Header_h */
