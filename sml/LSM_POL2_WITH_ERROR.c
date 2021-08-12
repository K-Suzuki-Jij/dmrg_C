#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void LSM_POL2_WITH_ERROR(double *Data_x,
                         double *Data_y,
                         double *Delta_y,
                         long data_num,
                         double *coeff2,
                         double *coeff1,
                         double *coeff0,
                         double *delta_coeff2,
                         double *delta_coeff1,
                         double *delta_coeff0,
                         double *R2
                         ){
   
   if (data_num <= 3) {
      printf("Error in LSM_POL2_WITH_ERROR\n");
      printf("data_num=%ld\n", data_num);
      exit(1);
   }
   
   long i;
   double x0   = 0;
   double x1   = 0;
   double x2   = 0;
   double x3   = 0;
   double x4   = 0;
   double y1   = 0;
   double x1y1 = 0;
   double x2y1 = 0;
   double temp;
   
   for (i = 0; i < data_num; i++) {
      temp = 1.0/(Delta_y[i]*Delta_y[i]);
      x4 = x4 + Data_x[i]*Data_x[i]*Data_x[i]*Data_x[i]*temp;
      x3 = x3 + Data_x[i]*Data_x[i]*Data_x[i]*temp;
      x2 = x2 + Data_x[i]*Data_x[i]*temp;
      x1 = x1 + Data_x[i]*temp;
      x0 = x0 + temp;
      x2y1 = x2y1 + Data_x[i]*Data_x[i]*Data_y[i]*temp;
      x1y1 = x1y1 + Data_x[i]*Data_y[i]*temp;
      y1 = y1 + Data_y[i]*temp;
   }
   
   double det = -x2*x2*x2 + 2*x1*x2*x3 - x0*x3*x3 - x1*x1*x4 + x0*x2*x4;
   if (fabs(det) < pow(10,-15)) {
      printf("Error in LSM_POL2_WITH_ERROR\n");
      printf("det=%.15lf\n", det);
      exit(1);
   }
   
   double a11 = (-x1*x1+x0*x2)/det;
   double a12 = (x1*x2-x0*x3)/det;
   double a13 = (-x2*x2+x1*x3)/det;
   double a21 = (x1*x2-x0*x3)/det;
   double a22 = (-x2*x2+x0*x4)/det;
   double a23 = (x2*x3-x1*x4)/det;
   double a31 = (-x2*x2+x1*x3)/det;
   double a32 = (x2*x3-x1*x4)/det;
   double a33 = (-x3*x3+x2*x4)/det;
   
   double c0,c1,c2;
   
   c0 = a31*x2y1 + a32*x1y1 + a33*y1;
   c1 = a21*x2y1 + a22*x1y1 + a23*y1;
   c2 = a11*x2y1 + a12*x1y1 + a13*y1;
   
   *coeff0 = c0;
   *coeff1 = c1;
   *coeff2 = c2;
   
   //Calculate errors
   double disp_c0 = 0;
   double disp_c1 = 0;
   double disp_c2 = 0;
   double temp_c0;
   double temp_c1;
   double temp_c2;
   
   for (i = 0;i < data_num; i++) {
      temp = 1.0/Delta_y[i];
      temp_c0 = temp*(a31*Data_x[i]*Data_x[i] + a32*Data_x[i] + a33);
      temp_c1 = temp*(a21*Data_x[i]*Data_x[i] + a22*Data_x[i] + a23);
      temp_c2 = temp*(a11*Data_x[i]*Data_x[i] + a12*Data_x[i] + a13);
      disp_c0 = disp_c0 + temp_c0*temp_c0;
      disp_c1 = disp_c1 + temp_c1*temp_c1;
      disp_c2 = disp_c2 + temp_c2*temp_c2;
   }
   
   *delta_coeff0 = sqrt(disp_c0);
   *delta_coeff1 = sqrt(disp_c1);
   *delta_coeff2 = sqrt(disp_c2);

   double disp_f_y = 0;
   
   for (i = 0; i < data_num; i++) {
      temp = c2*Data_x[i]*Data_x[i] + c1*Data_x[i] + c0 - Data_y[i];
      disp_f_y = disp_f_y + temp*temp;
   }
   
   double mean_y = 0;
   
   for (i = 0; i < data_num; i++) {
      mean_y = mean_y + Data_y[i];
   }
   mean_y = mean_y/data_num;
   
   double disp_y = 0;
   
   for (i = 0; i < data_num; i++) {
      temp = mean_y - Data_y[i];
      disp_y = disp_y + temp*temp;
   }
   
   *R2 = 1 - disp_f_y/disp_y;

}
