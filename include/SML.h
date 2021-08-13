#ifndef SML_h
#define SML_h

typedef struct {
   
   long max_val;
   int max_row;
   int row_dim;
   int col_dim;
   
   int *Col;
   long *Row;
   double _Complex *Val;
   
} C_CRS1;

typedef struct {
   
   long max_val;
   int max_row;
   int row_dim;
   int col_dim;
   
   int *Col;
   long *Row;
   double *Val;
   
} CRS1;

typedef struct {
   
   long max_val;
   long max_row;
   long row_dim;
   long col_dim;
   
   long *Col;
   long *Row;
   double *Val;
   
} LCRS1;

typedef struct {
   
   long max_val;
   int max_col;
   int row_dim;
   int col_dim;
   
   long *Col;
   int *Row;
   double *Val;
   
} CCS1;

typedef struct {

   CRS1 *M;
   double acc;
   double *Vec;
   double *Out_Vec;
   int max_step;
   int p_threads;
   
} BOX_CG;

typedef struct {
   
   CRS1 *M;
   double ii_acc;
   double ii_diag_add;
   int ii_max_step;
   
   double cg_acc;
   int cg_max_step;

   double error;
   
   double *Eig_Vec;
   double *Eig_Val;
   
   int p_threads;
   char Type[100];

   
} BOX_II;

typedef struct {
   
   CRS1 *M;
   double ii_acc;
   double ii_diag_add;
   int ii_max_step;

   double cg_acc;
   int cg_max_step;
   
   double **Eig_Vec;
   double *Eig_Val;
   
   int eig_num;
   int p_threads;
   char Type[100];

} BOX_BLOCK_II;

typedef struct {

   CRS1 *M;

   double acc;
   double *Eig_Val;
   double **Eig_Vec;
   
   int max_step;
   int eig_num;
   int conv_eig_num;
   int p_threads;
   
} BOX_BLOCK_LOBPCG;

typedef struct {
   
   CRS1 *M;
   double *Eig_Vec;
   double *eig_val;
   double acc;
   int min_step;
   int max_step;
   int max_block_step;
   int p_threads;
   char Guess[100];
   char Type[100];

} BOX_LAN;

typedef struct {

   CRS1 *M;
   double *Eig_Vec;
   double *GS_Vec;
   double *eig_val;
   double acc;
   int min_step;
   int max_step;
   int p_threads;
   char Type[100];

   
} BOX_LAN_EX1;

typedef struct {

   CRS1 *M;
   int max_step;
   double acc;
   double *eig_val;
   double *Eig_Vec;
   int p_threads;

} BOX_LOBPCG;

long BINARY_SEARCH_INT1(int *Array, long imin, long imax, long target_val);
long BINARY_SEARCH_LINT1(long *Array, long imin, long imax, long target_val);
long BINOMIAL_COEFFICIENT(int n, int k);
void BLOCK_INVERSE_ITERATION(BOX_BLOCK_II *Box_II);
void BLOCK_LOBPCG(BOX_BLOCK_LOBPCG *Box);
void BUBBLE_SORT_INT1(int *Array, long dim);

void   CHECK_CRS1(CRS1 *M);
int    CHECK_SYMMETRY_CRS1(CRS1 *M, double zero, int p_threads);
void   CONJUGATE_GRADIENT(BOX_CG *Box);
void   CONJUGATE_GRADIENT_SYM(BOX_CG *Box);
double CONVERGE_CHECK(CRS1 *M, double *Eig_Vec, double eig_val, int p_threads);
double CONVERGE_CHECK_SYM(CRS1 *M, double *Eig_Vec, double eig_val, int p_threads);
void   COPY_CRS1(CRS1 *M, CRS1 *Copyed, int p_threads);
void   COPY_CCS1(CCS1 *M, CCS1 *Copyed, int p_threads);
void   COPY_DOUBLE1(double *V, double *Copyed, long dim, int p_threads);
void   COPY_INT1(int *V, int *Copyed, long dim, int p_threads);
void   COPY_LINT1(long *V, long *Copyed, long dim, int p_threads);
void   CRS_CRS_CCS_PRODUCT(CRS1 *M1, CRS1 *M2, CCS1 *M3, CRS1 *Out, CCS1 *Work, double *Work_Vec);

void DIAG_ADD_CRS1(CRS1 *M, double add, int p_threads);
int DELTA_FUNCTION(long i, long j);

void FREE_ARRAY_CHAR1(char *Matrix);
void FREE_ARRAY_CHAR2(char **Matrix,long row);
void FREE_ARRAY_CHAR3(char ***Matrix, long row, long col);
void FREE_ARRAY_DOUBLE1(double *Matrix);
void FREE_ARRAY_DOUBLE2(double **Matrix, long row);
void FREE_ARRAY_DOUBLE3(double ***Matrix, long row, long col);
void FREE_ARRAY_DOUBLE4(double ****Matrix, long row, long col, long col_2);
void FREE_ARRAY_DOUBLE5(double *****Matrix, long row, long col, long col_2, long col_3);
void FREE_ARRAY_DOUBLE6(double ******Matrix, long row, long col, long col_2, long col_3, long col_4);
void FREE_ARRAY_DOUBLE7(double *******Matrix, long row, long col, long col_2, long col_3, long col_4, long col_5);
void FREE_ARRAY_INT1(int *Matrix);
void FREE_ARRAY_INT2(int **Matrix, long row);
void FREE_ARRAY_INT3(int ***Matrix, long row, long col);
void FREE_ARRAY_INT4(int ****Matrix, long row, long col, long col_2);
void FREE_ARRAY_INT5(int *****Matrix, long row, long col, long col_2, long col_3);
void FREE_ARRAY_INT6(int ******Matrix, long row, long col, long col_2, long col_3, long col_4);
void FREE_ARRAY_LINT1(long *Matrix);
void FREE_ARRAY_LINT2(long **Matrix, long row);
void FREE_ARRAY_LINT3(long ***Matrix, long row, long col);
void FREE_ARRAY_LINT4(long ****Matrix, long row, long col, long col2);
void FREE_ARRAY_SINT1(short int *Matrix);
void FREE_ARRAY_SINT2(short int **Matrix, long row);
void FREE_ARRAY_SINT3(short int ***Matrix, long row, long col);
void FREE_ARRAY_SINT4(short int ****Matrix, long row, long col, long col_2);
void FREE_ARRAY_SINT5(short int *****Matrix, long row, long col, long col_2, long col_3);
void FREE_CCS1(CCS1 *Matrix);
void FREE_CCS2(CCS1 **Matrix, int row);
void FREE_CRS1(CRS1 *Matrix);
void FREE_LCRS1(LCRS1 *Matrix);
void FREE_CRS2(CRS1 **Matrix, int row);
void FREE_CRS3(CRS1 ***Matrix, int row1, int row2);
void FREE_CRS4(CRS1 ****Matrix, int row1, int row2, int row3);
int FIND_MIN_INT1(int *Array, int dim);
int FIND_MAX_INT1(int *Array, long dim);
int FIND_MAX_INT2(int **Array, long row1, long row2);
int FIND_MAX_INT3(int ***Array, long row1, long row2, long row3);
int FIND_MAX_INT4(int ****Array, long row1, long row2, long row3, long row4);
long FIND_MAX_LINT1(long *Array, long dim);
long FIND_MAX_LINT2(long **Array, long row1, long row2);
long FIND_MAX_LINT3(long ***Array, long row1, long row2, long row3);


char      *GET_ARRAY_CHAR1(long num);
char      **GET_ARRAY_CHAR2(long row, long col);
char      ***GET_ARRAY_CHAR3(long row, long col, long col_2);
double _Complex *GET_ARRAY_C_DOUBLE1(long num);
double    *GET_ARRAY_DOUBLE1(long num);
double    **GET_ARRAY_DOUBLE2(long row, long col);
double    ***GET_ARRAY_DOUBLE3(long row, long col, long col_2);
double    ****GET_ARRAY_DOUBLE4(long row, long col, long col_2, long col_3);
double    *****GET_ARRAY_DOUBLE5(long row, long col, long col_2, long col_3, long col_4);
double    ******GET_ARRAY_DOUBLE6(long row, long col, long col_2, long col_3, long col_4, long col_5);
double    *******GET_ARRAY_DOUBLE7(long row, long col, long col_2, long col_3, long col_4, long col_5, long col_6);
int       *GET_ARRAY_INT1(long num);
int       **GET_ARRAY_INT2(long row, long col);
int       ***GET_ARRAY_INT3(long row, long col, long col_2);
int       ****GET_ARRAY_INT4(long row, long col, long col_2, long col_3);
int       *****GET_ARRAY_INT5(long row, long col, long col_2, long col_3, long col_4);
int       ******GET_ARRAY_INT6(long row, long col, long col_2, long col_3, long col_4, long col_5);
long      *GET_ARRAY_LINT1(long num);
long      **GET_ARRAY_LINT2(long row, long col);
long      ***GET_ARRAY_LINT3(long row, long col, long col_2);
long      ****GET_ARRAY_LINT4(long row, long col, long col_2, long col_3);
short int *GET_ARRAY_SINT1(long num);
short int **GET_ARRAY_SINT2(long row, long col);
short int ***GET_ARRAY_SINT3(long row, long col, long col_2);
short int ****GET_ARRAY_SINT4(long row, long col, long col_2, long col_3);
short int *****GET_ARRAY_SINT5(long row, long col, long col_2, long col_3, long col_4);
C_CRS1    *GET_C_CRS1(int dim, long max);
CRS1      *GET_CRS1(int dim, long max);
CRS1      **GET_CRS2(int row, int dim, long max);
CRS1      ***GET_CRS3(int row1, int row2, int dim, long max);
CRS1      ****GET_CRS4(int row1, int row2, int row3, int dim, long max);
CCS1      *GET_CCS1(int dim, long max);
CCS1      **GET_CCS2(int row, int dim, long max);
CRS1      *GET_COPY_CRS1(CRS1 *M);
double    *GET_RAND_DOUBLE1(long num);

double INNER_PRODUCT(double *V1, double *V2, long dim, int p_threads);
void   INSERTION_SORT_LINT1_INT4_DOUBLE1(long *Int_Sort, int *Int_Array2, int *Int_Array3, int *Int_Array4, int *Int_Array5, double *D_Array2, long dim);
void   INVERSE_ITERATION(BOX_II *Box_II);
void   INVERSE_ITERATION_SYM(BOX_II *Box_II);
double L1_NORM(double *V, long dim, int p_threads);
double L2_NORM(double *V, long dim, int p_threads);
void   LANCZOS_SLOW_EX1(BOX_LAN_EX1 *Box);
void   LANCZOS_SLOW(BOX_LAN *Box);
void   LANCZOS_SLOW_SYM(BOX_LAN *Box);
void   LANCZOS(BOX_LAN *Box);
void   LANCZOS_SYM(BOX_LAN *Box);
void   LAPACK_DGEEV(double **Ham, int dim, double *Eigen_Value_Real, double *Eigen_Value_Img, double **Eigen_Vector_Real, double **Eigen_Vector_Img, int eigval_num, int eigvec_num);
void   LAPACK_DGESVD(double **M, int row, int col, double *Eigen_Value, double **Eigen_Vector_R, double **Eigen_Vector_L, int eigval_num, int eigvec_num_r, int eigvec_num_l);
void   LAPACK_DSYEV(double **Ham, int row, int col, double *Eigen_Value, double **Eigen_Vector, int eigval_num, int eigvec_num);
void   LAPACK_DSYEV_CRS1(CRS1 *Ham, double *Eigen_Value, double **Eigen_Vector, int eigval_num, int eigvec_num);
void   LSM_POL1_WITH_ERROR(double *Data_x, double *Data_y, double *Delta_y, long data_num, double *coeff1, double *coeff0, double *delta_coeff1, double *delta_coeff0, double *R2);
void   LSM_POL1(double *Data_x, double *Data_y, long data_num, double *coeff1, double *coeff0, double *delta_coeff1, double *delta_coeff0, double *R2);
void   LSM_POL2_WITH_ERROR(double *Data_x, double *Data_y, double *Delta_y, long data_num, double *coeff2, double *coeff1, double *coeff0, double *delta_coeff2, double *delta_coeff1, double *delta_coeff0, double *R2);
void   LSM_POL2(double *Data_x, double *Data_y, long data_num, double *coeff2, double *coeff1, double *coeff0, double *delta_coeff2, double *delta_coeff1, double *delta_coeff0, double *R2);

void MAKE_RAND_CRS1(CRS1 *M, long dim, long ele_num);
void MAKE_RAND_SYM_CRS1(CRS1 *M, long dim, long ele_num);
void MATRIX_CONSTAN_MULTIPLICATION_CRS1(CRS1 *M, double c, int p_threads);
void MATRIX_MATRIX_MATRIX_PRODUCT_CRS1(CRS1 *T_Matrix, CRS1 *Matrix, CRS1 *T_Dagger_Matrix, CRS1 *Out, CRS1 *Work);
void MATRIX_PRODUCT_CRS1(CRS1 *M1, CRS1 *M2, CRS1 *Out);
void MATRIX_SUM_CRS1(CRS1 *M1, CRS1 *M2, CRS1 *Out);
void MATRIX_TRANSPOSE_CRS1(CRS1 *M, CRS1 *Out);
void MATRIX_VECTOR_PRODUCT(CRS1 *M, double *V1, double *Out, int p_threads);
void MATRIX_VECTOR_PRODUCT_SYM(CRS1 *M, double *V1, double *Out, double **Temp, int p_threads);
void MINIMUM_RESIDUAL(BOX_CG *Box);
void MINIMUM_RESIDUAL_SYM(BOX_CG *Box);

void NORMALIZE(double *V, long dim, int p_threads);

void ORTHOGONALIZATION(double **Vector, long v_num, long dim, long start, int p_threads);

void PRINT_CRS1(CRS1 *Matrix, char Name[]);
void PRINT_CCS1(CCS1 *Matrix, char Name[]);
void PRINT_DOUBLE1(double *Vec, long dim, char Name[]);
void PRINT_INT1(int *Vec, long dim, char Name[]);
void PRINT_LINT1(long *Vec, long dim, char Name[]);
void PRINT_SINT1(short int *Vec, long dim, char Name[]);

void QUICK_SORT_LINT1(long *Target, long left, long right);
void QUICK_SORT_INT1_DOUBLE1(int *Target, double *D1, long left, long right);
void QUICK_SORT_INT1_SINT4(int *Target, short *A1, short *A2, short *A3, short *A4, int left, int right);
void QUICK_SORT_INT2_DOUBLE1(int *Target, int *A1, double *D2, long left, long right);
void QUICK_SORT_LINT2_DOUBLE1(long *Target, long *A1, double *D2, long left, long right);
void QUICK_SORT_STABLE_INT1_SINT4(int *Base, int *Target, short *A1, short *A2, short *A3, short *A4, int left, int right);
void QUICK_SORT_STABLE_INT2_SINT4(int *Base, int *Target, int *A_Int, short *A1, short *A2, short *A3, short *A4, int left, int right);
void QUICK_SORT_STABLE_INT3_SINT4(int *Base, int *Target, int *A1_Int, int *A2_Int, short *A1, short *A2, short *A3, short *A4, int left, int right);
void QUICK_SORT_STABLE_SINT7(int *Base, short *Target, short *A1, short *A2, short *A3, short *A4, short *A5, short *A6, int left, int right);

long RANDOM(long min, long max);
void READ_CRS1(CRS1 *M, char Name[]);

int SIGN(double a);
void SORT_COLUMN_CRS1(CRS1 *M, int p_threads);
void SWAP_DOUBLE(double *a, double *b);
void SWAP_INT(int *a, int *b);
void SWAP_LINT(long *a, long *b);
void SWAP_SINT(short *a, short *b);

double VECTOR_MATRIX_VECTOR_PRODUCT(double *V1, CRS1 *M, double *V2, int p_threads);

void WRITE_CRS1(CRS1 *M, char Name[]);

#endif /* SML_h */
