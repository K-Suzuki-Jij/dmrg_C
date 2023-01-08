//
//  Header.h
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/12.
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
   
   CRS1 **Ham;
   CRS1 **Sp_RE;
   CRS1 **Sm_RE;
   CRS1 **Sz_RE;
   
   CRS1 **Sp_LE;
   CRS1 **Sm_LE;
   CRS1 **Sz_LE;
   
   CRS1 *Sz_On;
   CRS1 *Sp_On;
   CRS1 *Sm_On;
   CRS1 *Ham_On;
   CRS1 *Sx_On;
   CRS1 *SzSz_On;
   CRS1 *SxSx_On;

   CRS1 **Sz;
   CRS1 **Sx;
   CRS1 **SzSz;
   CRS1 **SxSx;
   CRS1 **Sz_CF;
   CRS1 **Sx_CF;
   
   int **Tot_Sz;
   int *Dim;
   
} BLOCK;

typedef struct {
   
   CRS1 *Sp_On;
   CRS1 *Sm_On;
   CRS1 *Sz_On;
   CRS1 *Ham_On;
   
   CRS1 *Sp_RE_System;
   CRS1 *Sm_RE_System;
   CRS1 *Sz_RE_System;
   CRS1 *Ham_System;
   
   CRS1 *Sp_LE_System;
   CRS1 *Sm_LE_System;
   CRS1 *Sz_LE_System;
   
   CRS1 *Sp_RE_Enviro;
   CRS1 *Sm_RE_Enviro;
   CRS1 *Sz_RE_Enviro;
   CRS1 *Ham_Enviro;
   
   CRS1 *Sp_LE_Enviro;
   CRS1 *Sm_LE_Enviro;
   CRS1 *Sz_LE_Enviro;
   
   double SSD_LR;
   double SSD_RL;
   double SSD_LLLR;
   double SSD_LRRL;
   double SSD_RRRL;
   
   double J_xy;
   double J_z;
   double D_z;
   double h_z;
   
   char *BC;
   
} HAM_BOX;

//ALLOCATIONS
void FREE_BLOCK_MATRIX(BLOCK *Block, MODEL_1DXXZ_VF *Model);
void FREE_HAM_BOX(HAM_BOX **Box, MODEL_1DXXZ_VF *Model);
void FREE_WHOLE_BASIS_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis, MODEL_1DXXZ_VF *Model);
DMRG_BASIS *GET_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
BLOCK *GET_BLOCK_MATRIX(MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param);
DMRG_STATUS *GET_DMRG_STATUS(MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param);
HAM_BOX **GET_HAM_BOX(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, int LL_site, int RR_site);
CRS1 *GET_HAM_LLLR(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status, MODEL_1DXXZ_VF *Model);
CRS1 *GET_HAM_LRRL(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status, MODEL_1DXXZ_VF *Model);
CRS1 *GET_HAM_RRRL(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status, MODEL_1DXXZ_VF *Model);
CRS1 *GET_HAM_LLLRRRRL(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
DMRG_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void GET_GROUND_STATE(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);

//EXPECTATIONS
void EXPECTATION_INTERSITE_Q2(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, int tot_sz, double **Temp_V1, double **Temp_V2, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status);
void EXPECTATION_VALUES(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, double *Vec, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void TRANSFORM_MATRIX_FOR_EXPECTATION_VALUES(BLOCK *Block_System, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);

//MAIN
void DIAGONALIZE_SUPERBLOCK(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void DMRG(MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param);
double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DXXZ_VF *Model);
int FIND_SITE_SZ(int basis, int spin);
void MAKE_ELEMENT_HAM_LLLR(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box);
void MAKE_ELEMENT_HAM_LLLRRRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int *Inv, HAM_BOX *Box);
void MAKE_ELEMENT_HAM_LRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box);
void MAKE_ELEMENT_HAM_RRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box);
void MAKE_ONSITE_HAM(MODEL_1DXXZ_VF *Model, CRS1 *M);
void MAKE_OP_EDGE(BLOCK *Block, MODEL_1DXXZ_VF *Model);
void RENORMALIZE(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);
void TRANSFORM_MATRIX(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DXXZ_VF *Model, DMRG_TIME *Dmrg_Time);
//STATUS
void OUTPUT_ENERGY(MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_INTERSITE_VALUES(double *Out, double *Onsite_Val, char Name[], MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_ONSITE_VALUES(double *Out, char Name[], MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status);
void OUTPUT_T_ERROR(DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DXXZ_VF *Model);
void PRINT_MEM_STATUS(MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param);
void PRINT_STATUS(DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time);

#endif /* Header_h */
