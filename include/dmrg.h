#ifndef dmrg_h
#define dmrg_h

#include "SML.h"

typedef struct {
   
   char Initial_Guess[1000];
   char Enviro_Copy[1000];
   char II_Type[1000];
   char Lan_Con[1000];
   int max_dim_system;
   int sweep;
   
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
   
   double sp_LL;
   
   int param_iter;
   int param_iter_now;
   
   int dim_LLLRRRRL_limit;
   
   
} DMRG_PARAMETER;

typedef struct {
   
   double diag;
   double inv_iter;
   double make_ham;
   double make_basis;
   double make_dens_mat;
   double trans_main;
   double trans_exp;
   double total;
   
} DMRG_TIME;


typedef struct {
   
   short *LL_LLLRRRRL;
   short *LR_LLLRRRRL;
   short *RL_LLLRRRRL;
   short *RR_LLLRRRRL;
   int *Inv_LLLRRRRL;
   int dim_LLLRRRRL;

   int dim_LL;
   int dim_RR;
   int dim_onsite;
   
   //Model dependent
   int tot_sz_LLLRRRRL;
   int tot_ele_LLLRRRRL;
   int tot_ele_1_LLLRRRRL;
   int tot_ele_2_LLLRRRRL;
   int tot_parity_LLLRRRRL;
   
   short *LL_LLLR;
   short *LR_LLLR;
   int **Inv_LLLR;
   int dim_LLLR;
   
   short *LR_LRRL;
   short *RL_LRRL;
   int **Inv_LRRL;
   int dim_LRRL;
   
   short *RR_RRRL;
   short *RL_RRRL;
   int **Inv_RRRL;
   int dim_RRRL;
   
   int *Sum_LLLR;
   
   //Model dependent
   int *Tot_Sz_LLLR;
   int *Tot_Ele_LLLR;
   int *Tot_Ele_1_LLLR;
   int *Tot_Ele_2_LLLR;
   int *Tot_Parity_LLLR;
   
   
} DMRG_BASIS;

typedef struct {
   
   short *LL_LLLR;
   short *LR_LLLR;
   int **Inv_LLLR;
   int dim_LLLR;
   
} DMRG_BASIS_LLLR;

typedef struct {
   
   short *LR_LRRL;
   short *RL_LRRL;
   int **Inv_LRRL;
   int dim_LRRL;
   
} DMRG_BASIS_LRRL;

typedef struct {
   
   short *RR_RRRL;
   short *RL_RRRL;
   int **Inv_RRRL;
   int dim_RRRL;
   
} DMRG_BASIS_RRRL;

typedef struct {
   
   short *LL_LLLRRRRL;
   short *LR_LLLRRRRL;
   short *RL_LLLRRRRL;
   short *RR_LLLRRRRL;
   int *Inv_LLLRRRRL;
   double *Val_LLLRRRRL;
   int elem_num;
   
} DMRG_A_BASIS;

typedef struct {
   
   short **LL_LLLRRRRL;
   short **LR_LLLRRRRL;
   short **RL_LLLRRRRL;
   short **RR_LLLRRRRL;
   int *Inv_LLLRRRRL;
   int *Dim;
   int max_dim;
   int dim_RR;
   int dim_onsite;
   
} DMRG_WHOLE_BASIS_Q1;


typedef struct {
   
   short ***LL_LLLRRRRL;
   short ***LR_LLLRRRRL;
   short ***RL_LLLRRRRL;
   short ***RR_LLLRRRRL;
   int *Inv_LLLRRRRL;
   int **Dim;
   int max_dim;
   int dim_RR;
   int dim_onsite;
   
} DMRG_WHOLE_BASIS_Q2;

typedef struct {
   
   short ****LL_LLLRRRRL;
   short ****LR_LLLRRRRL;
   short ****RL_LLLRRRRL;
   short ****RR_LLLRRRRL;
   int *Inv_LLLRRRRL;
   int ***Dim;
   int max_dim;
   int dim_RR;
   int dim_onsite;
   
} DMRG_WHOLE_BASIS_Q3;

typedef struct {
   
   short *****LL_LLLRRRRL;
   short *****LR_LLLRRRRL;
   short *****RL_LLLRRRRL;
   short *****RR_LLLRRRRL;
   int *Inv_LLLRRRRL;
   int ****Dim;
   int max_dim;
   int dim_RR;
   int dim_onsite;
   
} DMRG_WHOLE_BASIS_Q4;


typedef struct {
   
   int LL;
   int LR;
   int RR;
   int RL;
   int row;
   int dim_RR;
   int dim_onsite;
   
} DMRG_BASIS_ONSITE;

typedef struct {
   
   int *Dim_Block;
   int *Q_Number1_Block;
   int *Q_Number2_Block;
   int *Q_Number3_Block;
   int block_num;
   
} DMRG_DENSITY_MATRIX_INFO;

typedef struct {
   
   CCS1 *Trans_Matrix;
   CRS1 *Trans_Matrix_Dagger;
   double tr_error;
   double *Val_DM_Dist;
   double sum_val_dm;
   int dim_LLLR;
   int dim_renorm;
   
   int *Q_Number1;
   int *Q_Number2;
   int *Q_Number3;
   
} DMRG_SYSTEM_INFO;

typedef struct {
   
   char BC[1000];
   char Enviro_Copy[100];
   int LL_site;
   int RR_site;
   
   int max_dim_system;
   int dim_LLLRRRRL;
   int dim_LL;
   int dim_onsite;
   int dim_RR;
   int spin;
   int tot_sz;
   int tot_ele;
   int tot_ele_1;
   int tot_ele_2;
   int tot_parity;
   double *GS_Vec;
   double gs_val;
   double gs_error;
   
   double percent_LL;
   
   int sweep_now;
   int param_iter_now;
   int tot_iter_now;
   int sweep;
   int param_iter;
   int tot_iter;
   
} DMRG_STATUS;

typedef struct {
   
   double *Eig_Vec;
   double *eig_val;
   double acc;
   int min_step;
   int max_step;
   int p_threads;
   char Guess[100];
   char Type[100];
   CRS1 *M_LLLR;
   CRS1 *M_LRRL;
   CRS1 *M_LRRL_Sign;
   CRS1 *M_RRRL;
   int *Ele_RR;
   DMRG_BASIS *Dmrg_Basis;
   
} DMRG_BOX_LAN;

typedef struct {

   double acc;
   double diag_val;
   double *Vec;
   double *Out_Vec;
   int max_step;
   int p_threads;
   CRS1 *M_LLLR     ;
   CRS1 *M_LRRL     ;
   CRS1 *M_LRRL_Sign;
   CRS1 *M_RRRL     ;
   int *Ele_RR;
   DMRG_BASIS *Dmrg_Basis;
   
} DMRG_BOX_CG;

typedef struct {
   
   double ii_acc;
   double ii_diag_add;
   int ii_max_step;
   
   double cg_acc;
   int cg_max_step;
   
   double *Eig_Vec;
   double *eig_val;
   
   double error;
   
   int p_threads;
   
   CRS1 *M_LLLR;
   CRS1 *M_LRRL;
   CRS1 *M_LRRL_Sign;
   CRS1 *M_RRRL;
   int *Ele_RR;
   DMRG_BASIS *Dmrg_Basis;

} DMRG_BOX_II;

void DMRG_MAKE_ELEM_LL_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_LR_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_RR_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_RL_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_LLLR_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_LL, CRS1 *M_LR, double coeef, int *LL_Ele,                           int *elem_num, char Type[], char Sign_Flag[]);
void DMRG_MAKE_ELEM_LLRL_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_LL, CRS1 *M_RL, double coeef, int *Ele_LL, int *Ele_RR, int *Ele_LR, int *elem_num, char Type[], char Sign_Flag[]);
void DMRG_MAKE_ELEM_LLRR_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_LL, CRS1 *M_RR, double coeef, int *Ele_LL, int *Ele_LR,              int *elem_num, char Type[], char Sign_Flag[]);
void DMRG_MAKE_ELEM_LRRL_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_LR, CRS1 *M_RL, double coeef, int *Ele_RR, int *Ele_LR,              int *elem_num, char Type[], char Sign_Flag[]);
void DMRG_MAKE_ELEM_LRRR_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_LR, CRS1 *M_RR, double coeef, int *Ele_LR,                           int *elem_num, char Type[], char Sign_Flag[]);
void DMRG_MAKE_ELEM_RRRL_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_RR, CRS1 *M_RL, double coeef, int *Ele_RR,                           int *elem_num, char Type[], char Sign_Flag[]);
void DMRG_MAKE_ELEM_ZERO_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, int *elem_num);

void DMRG_MAKE_ELEM_LL_LLLR(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LLLR,  CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_LR_LLLR(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LLLR, CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_LLLR_LLLR(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LLLR, CRS1 *M_LL, CRS1 *M_LR, double coeef, int *LL_Ele, int *elem_num, char Type[], char Sign_Flag[]);

void DMRG_MAKE_ELEM_LRRL_LRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LRRL, CRS1 *M_LR, CRS1 *M_RL, double coeef, int *Ele_LR, int *elem_num, char Type[], char Sign_Flag[]);

void DMRG_MAKE_ELEM_RL_RRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_RRRL, CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_RR_RRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_RRRL, CRS1 *M, double coeef, int *elem_num);
void DMRG_MAKE_ELEM_RRRL_RRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_RRRL, CRS1 *M_RR, CRS1 *M_RL, double coeef, int *Ele_RR, int *elem_num, char Type[], char Sign_Flag[]);

void DMRG_FREE_A_BASIS(DMRG_A_BASIS **A_Basis, int p_threads);
void DMRG_DIAGONALIZE_SUPERBLOCK(CRS1 *Ham, DMRG_TIME *Time, DMRG_PARAMETER *Param, DMRG_STATUS *Dmrg_Status, int p_threads);
void DMRG_FREE_DENSITY_MATRIX_INFO(DMRG_DENSITY_MATRIX_INFO *Info);
void DMRG_MAKE_LLLR_OP_LLLR(CRS1 *M_LL, CRS1 *M_LR, double coeef, int *Ele_LL, char Sign[], char Type[], DMRG_BASIS_LLLR *Dmrg_Basis, CRS1 *Out);
void DMRG_MAKE_LR_OP_LLLR(CRS1 *M_On, int *Ele_LL, char Sign[], DMRG_BASIS_LLLR *Dmrg_Basis, CRS1 *Out);
void DMRG_FREE_SYSTEM_INFO(DMRG_SYSTEM_INFO *Info);
void DMRG_MAKE_LL_OP_LLLR(CRS1 *M_LL, DMRG_BASIS_LLLR *Dmrg_Basis, CRS1 *Out);
void DMRG_EXTEND_AND_REDUCE(CRS1 *M_LL, CRS1 *M_On, CRS1 *Out, int *Ele_LL, char Type[], char Sign[], CCS1 *T_M, CRS1 *T_MD, CRS1 *Work_CRS, CCS1 *Work_CCS, double *Work_Vec, DMRG_BASIS_LLLR *Dmrg_Basis);
void DMRG_EXPECTATION_ONSITE(CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, double *Vec, double *Temp_V, int p_threads, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void DMRG_EXPECTATION_INTERSITE_Q0(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, double *Temp_V1, double *Temp_V2, int p_threads, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status);
void DMRG_FREE_BASIS_LLLRRRRL(DMRG_BASIS *Basis);
void DMRG_FREE_BASIS_LLLR(DMRG_BASIS *Basis, int dim_LL);

void DMRG_V_M_LL_Q0(CRS1 *M_LL, double *V, double *Out_V, int dim, int p_threads, DMRG_BASIS *Dmrg_Basis);
void DMRG_V_M_LR_Q0(CRS1 *M_LR, double *V, double *Out_V, int dim, int p_threads, DMRG_BASIS *Dmrg_Basis);
void DMRG_V_M_RR_Q0(CRS1 *M_RR, double *V, double *Out_V, int dim, int p_threads, DMRG_BASIS *Dmrg_Basis);
void DMRG_V_M_RL_Q0(CRS1 *M_RL, double *V, double *Out_V, int dim, int p_threads, DMRG_BASIS *Dmrg_Basis);
void DMRG_V_M_LL_Q1(CRS1 *M_LL, int qn_out, double *V, double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q1 *Dmrg_W_Basis);
void DMRG_V_M_LR_Q1(CRS1 *M_LR, int qn_out, double *V, int *Ele_LL, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q1 *Dmrg_W_Basis);
void DMRG_V_M_RL_Q1(CRS1 *M_RL, int qn_out, double *V, int *Ele_LL, int *Ele_LR, int *Ele_RR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q1 *Dmrg_W_Basis);
void DMRG_V_M_RR_Q1(CRS1 *M_RR, int qn_out, double *V, int *Ele_LL, int *Ele_LR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q1 *Dmrg_W_Basis);
void DMRG_V_M_LL_Q2(CRS1 *M_LL, int qn1_out, int qn2_out, double *V, double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis);
void DMRG_V_M_LR_Q2(CRS1 *M_LR, int qn1_out, int qn2_out, double *V, int *Ele_LL, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis);
void DMRG_V_M_RL_Q2(CRS1 *M_RL, int qn1_out, int qn2_out, double *V, int *Ele_LL, int *Ele_LR, int *Ele_RR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis);
void DMRG_V_M_RR_Q2(CRS1 *M_RR, int qn1_out, int qn2_out, double *V, int *Ele_LL, int *Ele_LR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis);
void DMRG_V_M_LL_Q3(CRS1 *M_LL, int qn1_out, int qn2_out, int qn3_out, double *V, double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q3 *Dmrg_W_Basis);
void DMRG_V_M_LR_Q3(CRS1 *M_LR, int qn1_out, int qn2_out, int qn3_out, double *V, int *Ele_LL, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q3 *Dmrg_W_Basis);
void DMRG_V_M_RL_Q3(CRS1 *M_RL, int qn1_out, int qn2_out, int qn3_out, double *V, int *Ele_LL, int *Ele_LR, int *Ele_RR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q3 *Dmrg_W_Basis);
void DMRG_V_M_RR_Q3(CRS1 *M_RR, int qn1_out, int qn2_out, int qn3_out, double *V, int *Ele_LL, int *Ele_LR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q3 *Dmrg_W_Basis);
void DMRG_V_M_RL_Q4(CRS1 *M_RL, int qn1_out, int qn2_out, int qn3_out, int qn4_out, double *V, int *Ele_LL, int *Ele_LR, int *Ele_RR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q4 *Dmrg_W_Basis);
void DMRG_V_M_RR_Q4(CRS1 *M_RR, int qn1_out, int qn2_out, int qn3_out, int qn4_out, double *V, int *Ele_LL, int *Ele_LR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q4 *Dmrg_W_Basis);
double **DMRG_GET_DENSITY_MATRIX(int block, int *Dim_Block, DMRG_BASIS *Dmrg_Basis, double *Vec, int p_threads);
double DMRG_SSD_COEFF(int LL_site, int tot_site, char BC[], char Block_Name[], char Inter_Name[]);

void DMRG_GET_BASIS_LLLRRRRL(long dim_LLLRRRRL, int dim_onsite, int dim_LL, int dim_RR, DMRG_BASIS *Basis);
void DMRG_GET_BASIS_LLLR(int dim_onsite, int dim_LL, DMRG_BASIS *Basis);
void DMRG_GET_BASIS_RRRL(int dim_RR, int dim_onsite, DMRG_BASIS *Basis);
void DMRG_GET_BASIS_LRRL(int dim_onsite, DMRG_BASIS *Basis);
DMRG_SYSTEM_INFO *DMRG_GET_SYSTEM_INFO_Q1(int *Q_Number, DMRG_BASIS *Dmrg_Basis, double *Vec, int max_dim, int p_threads, DMRG_TIME *Dmrg_Time);
DMRG_SYSTEM_INFO *DMRG_GET_SYSTEM_INFO_Q2(int *Q_Number1, int *Q_Number2, DMRG_BASIS *Dmrg_Basis, double *Vec, int max_dim, int p_threads, DMRG_TIME *Dmrg_Time);
DMRG_SYSTEM_INFO *DMRG_GET_SYSTEM_INFO_Q3(int *Q_Number1, int *Q_Number2, int *Q_Number3, DMRG_BASIS *Dmrg_Basis, double *Vec, int max_dim, int p_threads, DMRG_TIME *Dmrg_Time);
DMRG_DENSITY_MATRIX_INFO *DMRG_GET_DENSITY_MATRIX_INFO_Q1(int *Q_Number_LLLR, int dim_LLLR);
DMRG_DENSITY_MATRIX_INFO *DMRG_GET_DENSITY_MATRIX_INFO_Q2(int *Q_Number1_LLLR, int *Q_Number2_LLLR, int dim_LLLR);
DMRG_DENSITY_MATRIX_INFO *DMRG_GET_DENSITY_MATRIX_INFO_Q3(int *Q_Number1_LLLR, int *Q_Number2_LLLR, int *Q_Number3_LLLR, int dim_LLLR);
DMRG_A_BASIS **DMRG_GET_A_BASIS(int max_elem_num, int p_threads);
void DMRG_FREE_MEMORY(DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System);
void DMRG_QUICK_SORT_BASIS_Q3(short *Array_1, short *Array_2, short *Array_3, short *Array_4, short *Array_5, short *Array_6, short *Array_7, int left, int right);
void DMRG_QUICK_SORT_BASIS_Q2(short *Array_1, short *Array_2, short *Array_3, short *Array_4, short *Array_5, short *Array_6, int left, int right);
void DMRG_QUICK_SORT_BASIS_Q1(short *Array_1, short *Array_2, short *Array_3, short *Array_4, short *Array_5, int left, int right);
void DMRG_MATRIX_VECTOR_PRODUCT_OBC(CRS1 *M_LLLR, CRS1 *M_LRRL, CRS1 *M_LRRL_Sign, CRS1 *M_RRRL, double *V, double *Out_V, int *Ele_RR, int p_threads, DMRG_BASIS *Dmrg_Basis);
void DMRG_LANCZOS_SLOW_OBC(DMRG_BOX_LAN *Box);
void DMRG_CONJUGATE_GRADIENT_OBC(DMRG_BOX_CG *Box);
void DMRG_INVERSE_ITERATION_OBC(DMRG_BOX_II *Box_II);
double DMRG_CONVERGE_CHECK_OBC(double eig_val, CRS1 *M_LLLR, CRS1 *M_LRRL, CRS1 *M_LRRL_Sign, CRS1 *M_RRRL, double *Eig_V, double *Temp_V, int *Ele_RR, int p_threads, DMRG_BASIS *Dmrg_Basis);
void DMRG_DIAGONALIZE_SUPERBLOCK_SYM(CRS1 *Ham, DMRG_TIME *Time, DMRG_PARAMETER *Param, DMRG_STATUS *Dmrg_Status, int p_threads);
void DMRG_RE_ALLOCATE_INV_LLLRRRRL(DMRG_BASIS *Basis, int dim_LL, int dim_RR, int dim_onsite, int p_threads);

void DMRG_TRANS_MAT_ONE(CRS1 *M_On, CRS1 **Out, int *Dim, int *Dim_LLLR, CCS1 **TM, CRS1 **TM_D, short **LL_LLLR, short **LR_LLLR, int ***Inv_LLLR, int tot_site, int p_threads);
void DMRG_TRANS_MAT_TWO(CRS1 *M_On, CRS1 **Out, int *Dim_LL, int *Dim_LLLR, CCS1 **TM, CRS1 **TM_D, short **LL_LLLR, short **LR_LLLR, int ***Inv_LLLR, int cf_origin, int tot_site, int p_threads);
void DMRG_TRANS_MAT_C_C(CRS1 *M_On_1, CRS1 *M_On_2, CRS1 **Out, int *Dim_LL, int *Dim_LLLR, CCS1 **TM, CRS1 **TM_D, short **LL_LLLR, short **LR_LLLR, int ***Inv_LLLR, int **Tot_Ele_LL, int tot_site);
#endif /* dmrg_h */

