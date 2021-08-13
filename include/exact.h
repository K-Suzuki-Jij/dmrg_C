//
//  exact.h
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/09.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#ifndef exact_h
#define exact_h

#include "SML.h"

typedef struct {
   
   int diag_num;
   int est_max_row_elem_num;
   
   //Lanczos or LOBPCG
   char Diag_Method[1000];
   double diag_acc;
   int diag_min_step;
   int diag_max_step;
   
   //Conjugate gradient Method or Minimal Residual
   double cg_acc;
   int cg_max_step;
   
   //Inverse iteration
   double inv_iter_acc;
   double inv_iter_diag_add;
   int inv_iter_max_step;
   
   int param_iter;
   int param_iter_now;
   
} EXACT_PARAMETER;

typedef struct {
   
   double diag;
   double inv_iter;
   double make_ham;
   double find_dim;
   double make_basis;
   double exp_values;
   double cf;
   double sc;
   double total;
   
} EXACT_TIME;


typedef struct {
   
   int dim;
   int diag_num;
   double **Vector;
   double *Value;
   double *Error;
   double mem_ham;
   CRS1 *Ham;
   
} EXACT_HAM_INFO;

typedef struct {
   
   int dim;
   long *Basis;
   int tot_sz;
   int tot_ele;
   int tot_parity;
   
} EXACT_BASIS_INFO;

typedef struct {
   
   long *Basis;
   long elem_num;
   long *Check;
   double *Val;
   
   int max_row;
   int dim_onsite;
   long dim_whole;
   
} EXACT_A_BASIS;

typedef struct {
   
   int *Dim;
   long **Basis;
   int max_dim;
   
} EXACT_WHOLE_BASIS_Q1;

typedef struct {
   
   int **Dim;
   long ***Basis;
   int max_dim;
   
} EXACT_WHOLE_BASIS_Q2;

typedef struct {
   
   int ***Dim;
   long ****Basis;
   int max_dim;
   
} EXACT_WHOLE_BASIS_Q3;

EXACT_A_BASIS **EXACT_GET_A_BASIS(int p_threads, int max_row);
int EXACT_FIND_SITE_STATE(long basis, int site, int dim_onsite);
void EXACT_MAKE_ELEM_INTER(long basis, int site1, int site2, int dim_onsite, CRS1 *M1, CRS1 *M2, long *elem_num, double coeef, int sign, EXACT_A_BASIS *A_Basis);
void EXACT_MAKE_ELEM_ON(long basis, int site, int dim_onsite, CRS1 *M_On, long *elem_num, double coeef, EXACT_A_BASIS *A_Basis);
void EXACT_FREE_A_BASIS(EXACT_A_BASIS **A_Basis, int p_threads);
void EXACT_DIAGONALIZE_HAMILTONIAN(EXACT_HAM_INFO *Ham_Info, EXACT_PARAMETER *Param, EXACT_TIME *Time, int p_threads);
void EXACT_V_M_Q0(CRS1 *M_On, double *Vec, double *Out_Vec, int dim_onsite, int site, int p_threads, EXACT_BASIS_INFO *Basis_Info);
void EXACT_EXPECTATION_ONSITE(CRS1 *M_On, double *Out, double *Vec, double *T_Vec, int dim_onsite, int tot_site, int p_threads, EXACT_BASIS_INFO *Basis_Info);
void EXACT_EXPECTATION_INTERSITE_Q0(CRS1 *M_O, CRS1 *M_R, int start, int end, double *Out, double *Vec, double *T_Vec1, double *T_Vec2, int dim_onsite, int p_threads, EXACT_BASIS_INFO *Basis_Info);
void EXACT_V_M_Q1(CRS1 *M_On, int qn_out, double *Vec, int qn_in, double *Out_Vec, int dim_onsite, int site, int p_threads, EXACT_WHOLE_BASIS_Q1 *W_Basis);
void EXACT_V_M_Q2(CRS1 *M_On, int qn1_out, int qn2_out, double *Vec, int qn1_in, int qn2_in, double *Out_Vec, char Sign_Flag[], int *N_Ele, int dim_onsite, int op_site, int p_threads, EXACT_WHOLE_BASIS_Q2 *W_Basis);
void EXACT_V_M_Q3(CRS1 *M_On, int qn1_out, int qn2_out, int qn3_out, double *Vec, int qn1_in, int qn2_in, int qn3_in, double *Out_Vec, char Sign_Flag[], int *N_Ele, int dim_onsite, int op_site, int p_threads, EXACT_WHOLE_BASIS_Q3 *W_Basis);
double EXACT_SSD_COEFF(int site, int tot_site, char BC[], char Inter_Name[]);
#endif /* exact_h */
