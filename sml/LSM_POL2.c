#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void LSM_POL2(double *Data_x,
              double *Data_y,
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
      printf("Error in LSM_POL2\n");
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
      x4 = x4 + Data_x[i]*Data_x[i]*Data_x[i]*Data_x[i];
      x3 = x3 + Data_x[i]*Data_x[i]*Data_x[i];
      x2 = x2 + Data_x[i]*Data_x[i];
      x1 = x1 + Data_x[i];
      x2y1 = x2y1 + Data_x[i]*Data_x[i]*Data_y[i];
      x1y1 = x1y1 + Data_x[i]*Data_y[i];
      y1 = y1 + Data_y[i];
   }
   x0 = data_num;
   
   double det = -x2*x2*x2 + 2*x1*x2*x3 - x0*x3*x3 - x1*x1*x4 + x0*x2*x4;
   if (fabs(det) < pow(10,-15)) {
      printf("Error in LSM_POL2\n");
      printf("det=%.15lf\n",det);
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
   double delta_y = 0;
   double disp_f_y = 0;
   for (i = 0; i < data_num; i++) {
      temp = c2*Data_x[i]*Data_x[i] + c1*Data_x[i] + c0 - Data_y[i];
      delta_y = delta_y + temp*temp;
   }
   
   disp_f_y = delta_y;
   delta_y = delta_y/(data_num - 3);
   *delta_coeff0 = sqrt(delta_y*(a31*a31*x4 + 2*a31*a32*x3 + (a32*a32 + 2*a31*a33)*x2 + 2*a32*a33*x1 + a33*a33*x0));
   *delta_coeff1 = sqrt(delta_y*(a21*a21*x4 + 2*a21*a22*x3 + (a22*a22 + 2*a21*a23)*x2 + 2*a22*a23*x1 + a23*a23*x0));
   *delta_coeff2 = sqrt(delta_y*(a11*a11*x4 + 2*a11*a12*x3 + (a12*a12 + 2*a11*a13)*x2 + 2*a12*a13*x1 + a13*a13*x0));
   
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
