#ifndef model_h
#define model_h

#include "SML.h"

typedef struct {
   
   char BC[1000];
   int p_threads;
   int spin;
   int tot_site;
   int tot_sz;
   int dim_onsite;
   int cf_origin;
   
   double J_xy;
   double J_z ;
   double D_z ;
   double h_z ;
   
} MODEL_1DXXZ_VF;

typedef struct {
   
   char BC[1000];
   int p_threads;
   int spin_loc;
   int tot_site;
   int tot_sz;
   int tot_ele;
   int dim_onsite;
   int dim_lspin;
   int dim_charge;
   int cf_origin;
   double avg_range;
   
   double t;
   double J;
   double D_z;
   double I_xy;
   double I_z;
   double h_z;
   double mu;
   
} MODEL_1DKLM_VF;

typedef struct {
   
   char BC[1000];
   int p_threads;
   int spin_loc;
   int tot_site;
   int tot_sz;
   int tot_ele_1;
   int tot_ele_2;
   int dim_onsite;
   int dim_lspin;
   int dim_charge;
   int dim_ccsl_onsite;
   int cf_origin;
   double avg_range;
   
   double t;
   double J;
   double D_z;
   double I_xy;
   double I_z;
   double h_z;
   double mu;
   
} MODEL_1DTKLM_VF;

typedef struct {
   
   char BC[1000];
   int p_threads;
   int spin_loc;
   int tot_site;
   int tot_parity;
   int tot_ele;
   int dim_onsite;
   int dim_lspin;
   int dim_charge;
   int cf_origin;
   double avg_range;
   
   double t;
   double J;
   double D_z;
   double I_xy;
   double I_z;
   double h_xc;
   double h_xl;
   double mu;
   
} MODEL_1DKLM_TVF;

typedef struct {
   
   char BC[1000];
   int p_threads;
   int tot_site;
   int tot_sz;
   int tot_ele;
   int dim_onsite;
   int cf_origin;
   double avg_range;
   
   double t1;
   double t2;
   double U;
   double V;
   double h_z;
   double mu;
   
   
} MODEL_1DHUBBARD_VF;

typedef struct {
   
   int max_dim;
   int dim_cc_1;
   int dim_cc_2;
   int dim_c_c;
   int dim_tot;
   
   int *CC_1_Num;
   int *CC_2_Num;
   int *C_Num1;
   int *C_Num2;
   int *CC_Parity;
   int *C_Parity;
   char **Row_Name;
   char **Col_Name;
   double ***Mat;
   double *F_Norm;
   
   //For SC correlations
   CRS1 **CC_Onsite   ;
   CRS1 **C_Onsite    ;
   CRS1 ***CC_Enviro  ;
   CRS1 ****C_C_Enviro;
   
} SC_MAT_1DKLM_TVF;

typedef struct {
   
   int max_dim;
   int dim_cc_1;
   int dim_cc_2;
   int dim_c_c;
   int dim_tot;
   
   int *CC_1_Num;
   int *CC_2_Num;
   int *C_Num1;
   int *C_Num2;
   int *CC_Sz;
   int *C_Sz;
   char **Row_Name;
   char **Col_Name;
   double ***Mat;
   double *F_Norm;
   
   //For SC correlations
   CRS1 **CC_Onsite   ;
   CRS1 **C_Onsite    ;
   CRS1 ***CC_Enviro  ;
   CRS1 ****C_C_Enviro;
   
} SC_MAT_1DKLM_VF;

typedef struct {
   
   int max_dim;
   int dim_cc_1;
   int dim_cc_2;
   int dim_c_c;
   int dim_tot;
   
   int *CC_1_Num;
   int *CC_2_Num;
   int *C_Num1;
   int *C_Num2;
   int *CC_Sz;
   int *C_Sz;
   char **Row_Name;
   char **Col_Name;
   double ***Mat;
   double *F_Norm;
   
   //For SC correlations
   CRS1 **CC_Onsite   ;
   CRS1 **C_Onsite    ;
   CRS1 ***CC_Enviro  ;
   CRS1 ****C_C_Enviro;
   
} SC_MAT_1DHUBBARD_VF;

typedef struct {
   
   int max_dim;
   int dim_ccsl;
   
   int *CCSL_Num;
   int *CCSL_Sz;
   int *CCSL_Ele_1;
   int *CCSL_Ele_2;
   char **Row_Name;
   char **Col_Name;
   double ***Mat;
   
} SC_MAT_1DTKLM_VF;

#endif /* model_h */
