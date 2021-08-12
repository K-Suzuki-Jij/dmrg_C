#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void LSM_POL1_WITH_ERROR(double *Data_x,
                         double *Data_y,
                         double *Delta_y,
                         long data_num,
                         double *coeff1,
                         double *coeff0,
                         double *delta_coeff1,
                         double *delta_coeff0,
                         double *R2
                         ){
   
   if (data_num <= 2) {
      printf("Error in LSM_POL1_WITH_ERROR\n");
      printf("data_num=%ld\n", data_num);
      exit(1);
   }

   long i;
   double temp;
   double x0,x1,x2,y1,x1y1;
   
   x2   = 0;
   x1   = 0;
   x0   = 0;
   x1y1 = 0;
   y1   = 0;

   for (i = 0; i < data_num; i++) {
      temp = 1.0/(Delta_y[i]*Delta_y[i]);
      x2 = x2 + Data_x[i]*Data_x[i]*temp;
      x1 = x1 + Data_x[i]*temp;
      x0 = x0 + temp;
      x1y1 = x1y1 + Data_x[i]*Data_y[i]*temp;
      y1 = y1 + Data_y[i]*temp;
   }
   
   double det = x2*x0 - x1*x1;
   if (fabs(det) < pow(10,-15)) {
      printf("Error in LSM_POL1_WITH_ERROR\n");
      printf("det=%.15lf\n",det);
      exit(1);
   }
   
   double a11 = x0/det;
   double a12 = -x1/det;
   double a21 = -x1/det;
   double a22 = x2/det;
   
   double c0,c1;
   
   c0 = a21*x1y1 + a22*y1;
   c1 = a11*x1y1 + a12*y1;
   
   *coeff0 = c0;
   *coeff1 = c1;
   
   //Calculate errors
   double disp_c0 = 0;
   double disp_c1 = 0;
   double temp_c0 = 0;
   double temp_c1 = 0;
   
   for (i = 0;i < data_num; i++) {
      temp = 1/Delta_y[i];
      temp_c0 = temp*(a21*Data_x[i] + a22);
      temp_c1 = temp*(a11*Data_x[i] + a12);
      disp_c0 = disp_c0 + temp_c0*temp_c0;
      disp_c1 = disp_c1 + temp_c1*temp_c1;
   }
   
   *delta_coeff0 = sqrt(disp_c0);
   *delta_coeff1 = sqrt(disp_c1);
   
   double disp_f_y = 0;
   for (i = 0; i < data_num; i++) {
      temp     = c1*Data_x[i] + c0 - Data_y[i];
      disp_f_y = disp_f_y + temp*temp;
   }
   
   double mean_y = 0;
   for (i = 0; i < data_num; i++) {
      mean_y = mean_y + Data_y[i];
   }
   mean_y = mean_y/data_num;
   
   double disp_y = 0;
   for (i = 0; i < data_num; i++) {
      temp   = mean_y - Data_y[i];
      disp_y = disp_y + temp*temp;
   }
   
   *R2 = 1 - disp_f_y/disp_y;
   
}
