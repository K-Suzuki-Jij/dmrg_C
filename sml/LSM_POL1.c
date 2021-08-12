//
//  LSM_POL1.c
//  
//
//  Created by 鈴木浩平 on 2018/05/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void LSM_POL1(double *Data_x,
              double *Data_y,
              long data_num,
              double *coeff1,
              double *coeff0,
              double *delta_coeff1,
              double *delta_coeff0,
              double *R2
              ){
   
   if (data_num <= 2) {
      printf("Error in LSM_POL1\n");
      printf("data_num=%ld\n", data_num);
      exit(1);
   }
   
   long i;
   double temp;
   double x0,x1,x2,y1,x1y1;
   
   x2   = 0;
   x1   = 0;
   x1y1 = 0;
   y1   = 0;
   
   for (i = 0; i < data_num; i++) {
      x2   = x2 + Data_x[i]*Data_x[i];
      x1   = x1 + Data_x[i];
      x1y1 = x1y1 + Data_x[i]*Data_y[i];
      y1   = y1 + Data_y[i];
   }
   x0 = data_num;
   
   double det = x2*x0 - x1*x1;
   if (fabs(det) < pow(10,-15)) {
      printf("Error in LSM_POL1\n");
      printf("det=%.15lf\n", det);
      exit(1);
   }
   
   double a11 =  x0/det;
   double a12 = -x1/det;
   double a21 = -x1/det;
   double a22 =  x2/det;
   
   double c0,c1;
   
   c0 = a21*x1y1 + a22*y1;
   c1 = a11*x1y1 + a12*y1;
   
   *coeff0 = c0;
   *coeff1 = c1;
   
   //Calculate errors
   double delta_y = 0;
   double disp_f_y = 0;
   
   for (i = 0; i < data_num; i++) {
      temp    = c1*Data_x[i] + c0 - Data_y[i];
      delta_y = delta_y + temp*temp;
   }
   disp_f_y = delta_y;
   delta_y = delta_y/(data_num - 2);
   
   *delta_coeff0 = sqrt(delta_y*(a21*a12*x2 + 2*a21*a22*x1 + a22*a22*x0));
   *delta_coeff1 = sqrt(delta_y*(a11*a11*x2 + 2*a11*a12*x1 + a12*a12*x0));
   
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
