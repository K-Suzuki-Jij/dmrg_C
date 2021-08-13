//
//  OUTPUT_SC_CORRELATIONS.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/30.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_SC_CORRELATIONS(SC_MAT_1DKLM_TVF *Sc_Mat, int site_start, int site_end, char File_Name[], char D1[], char D2[], MODEL_1DKLM_TVF *Model) {
   
   if (Sc_Mat->dim_tot <= 0) {
      return;
   }
   
   int i,j,k,count,r;
   int dim    = Sc_Mat->dim_tot;
   int length = site_end - site_start;
   double **Val_Real_Right = GET_ARRAY_DOUBLE2(length, dim);
   double **Val_Img_Right  = GET_ARRAY_DOUBLE2(length, dim);
   double **Val_Real_Left  = GET_ARRAY_DOUBLE2(length, dim);
   double **Val_Img_Left   = GET_ARRAY_DOUBLE2(length, dim);
   
   double ***Vec_Real_Right = GET_ARRAY_DOUBLE3(length, dim, dim);
   double ***Vec_Img_Right  = GET_ARRAY_DOUBLE3(length, dim, dim);
   double ***Vec_Real_Left  = GET_ARRAY_DOUBLE3(length, dim, dim);
   double ***Vec_Img_Left   = GET_ARRAY_DOUBLE3(length, dim, dim);
   double ***M_L            = GET_ARRAY_DOUBLE3(length, dim, dim);
   
   for (r = 0; r < length; r++) {
      for (i = 0; i < dim; i++) {
         for (j = 0; j < dim; j++) {
            M_L[r][i][j] = Sc_Mat->Mat[r][j][i];
         }
      }
   }
   
   for (r = 0; r < length; r++) {
      LAPACK_DGEEV(Sc_Mat->Mat[r], dim, Val_Real_Right[r], Val_Img_Right[r], Vec_Real_Right[r], Vec_Img_Right[r], dim, dim);
      LAPACK_DGEEV(M_L[r], dim, Val_Real_Left[r] , Val_Img_Left[r] , Vec_Real_Left[r] , Vec_Img_Left[r] , dim, dim);
   }
   
   FREE_ARRAY_DOUBLE3(M_L, length, dim);
   
   //Bubble Sort for Right Eigevalue
   
   for (r = 0; r < length; r++) {
      for (i = 0; i < dim - 1; i++) {
         count = 0;
         for (j = 0;j < dim - 1 - i; j++) {
            if ( fabs(Val_Real_Right[r][j]) < fabs(Val_Real_Right[r][j+1]) ) {
               
               SWAP_DOUBLE(&Val_Real_Right[r][j], &Val_Real_Right[r][j+1]);
               SWAP_DOUBLE(&Val_Img_Right[r][j] , &Val_Img_Right[r][j+1] );
               
               for (k = 0; k < dim; k++) {
                  SWAP_DOUBLE(&Vec_Real_Right[r][j][k], &Vec_Real_Right[r][j+1][k]);
                  SWAP_DOUBLE(&Vec_Img_Right[r][j][k] , &Vec_Img_Right[r][j+1][k] );
               }
               count = 1;
            }
         }
         if(count == 0){
            break;
         }
      }
   }
   
   //Bubble Sort for Left Eigenvalue
   for (r = 0; r < length; r++) {
      for (i = 0; i < dim - 1; i++) {
         count = 0;
         for (j = 0;j < dim - 1 - i; j++) {
            if ( fabs(Val_Real_Left[r][j]) < fabs(Val_Real_Left[r][j+1]) ) {
               
               SWAP_DOUBLE(&Val_Real_Left[r][j], &Val_Real_Left[r][j+1]);
               SWAP_DOUBLE(&Val_Img_Left[r][j] , &Val_Img_Left[r][j+1] );
               
               for (k = 0; k < dim; k++) {
                  SWAP_DOUBLE(&Vec_Real_Left[r][j][k], &Vec_Real_Left[r][j+1][k]);
                  SWAP_DOUBLE(&Vec_Img_Left[r][j][k] , &Vec_Img_Left[r][j+1][k] );
               }
               count = 1;
            }
         }
         if(count == 0){
            break;
         }
      }
   }
   
   char R_D1[300];
   char R_D2[300];
   char Val_D[300];
   char Vec_D[300];
   sprintf(R_D1 , "./result/%s"          , D1);
   sprintf(R_D2 , "./result/%s/%s"       , D1, D2);
   sprintf(Val_D, "./result/%s/%s/Value" , D1, D2);
   
   mkdir("./result", 0777);
   mkdir(R_D1, 0777);
   mkdir(R_D2, 0777);
   mkdir(Val_D, 0777);
   
   char File_Name_Val[300];
   char File_Name_Vec[300];
   char File_Basis_Row_Name[300];
   char File_Basis_Col_Name[300];
   FILE *file_val;
   FILE *file_vec;
   FILE *file_basis_row;
   FILE *file_basis_col;
   
   sprintf(File_Basis_Row_Name, "%s/%s_Row_Basis.txt", R_D2, File_Name);
   sprintf(File_Basis_Col_Name, "%s/%s_Col_Basis.txt", R_D2, File_Name);
   file_basis_row = fopen(File_Basis_Row_Name, "a+");
   file_basis_col = fopen(File_Basis_Col_Name, "a+");
   fprintf(file_basis_row ,"###N=%d,Ne=%d,LocSpin=%1.1lf,P=%d,BC=%s\n###t=%.1lf,J=%.3lf,I_xy=%.3lf,I_z=%.3lf,D_z=%.3lf,h_x=%.3lf,mu=%.5lf\n",
           Model->tot_site,
           Model->tot_ele,
           (double)Model->spin_loc/2.0,
           Model->tot_parity,
           Model->BC,
           Model->t,
           Model->J,
           Model->I_xy,
           Model->I_z,
           Model->D_z,
           Model->h_x,
           Model->mu
           );
   
   fprintf(file_basis_col ,"###N=%d,Ne=%d,LocSpin=%1.1lf,P=%d,BC=%s\n###t=%.1lf,J=%.3lf,I_xy=%.3lf,I_z=%.3lf,D_z=%.3lf,h_x=%.3lf,mu=%.5lf\n",
           Model->tot_site,
           Model->tot_ele,
           (double)Model->spin_loc/2.0,
           Model->tot_parity,
           Model->BC,
           Model->t,
           Model->J,
           Model->I_xy,
           Model->I_z,
           Model->D_z,
           Model->h_x,
           Model->mu
           );
   
   for (i = 0; i < dim; i++) {
      fprintf(file_basis_row, "%-3d  %-s\n", i, Sc_Mat->Row_Name[i]);
      fprintf(file_basis_col, "%-3d  %-s\n", i, Sc_Mat->Col_Name[i]);
   }
   fprintf(file_basis_row, "\n");
   fprintf(file_basis_col, "\n");
   fclose(file_basis_row);
   fclose(file_basis_col);
   
   
   int out_num;
   
   if (dim <= 5) {
      out_num = dim;
   }
   else {
      out_num = 5;
   }
   
   for (i = 0; i < out_num; i++) {
      sprintf(File_Name_Val, "%s/%s_Val_%d.txt", Val_D, File_Name, i);
      file_val = fopen(File_Name_Val, "a+");
      fprintf(file_val,"###N=%d,Ne=%d,LocSpin=%1.1lf,P=%d,BC=%s\n###t=%.1lf,J=%.3lf,I_xy=%.3lf,I_z=%.3lf,D_z=%.3lf,h_x=%.3lf,mu=%.5lf\n",
              Model->tot_site,
              Model->tot_ele,
              (double)Model->spin_loc/2.0,
              Model->tot_parity,
              Model->BC,
              Model->t,
              Model->J,
              Model->I_xy,
              Model->I_z,
              Model->D_z,
              Model->h_x,
              Model->mu
              );
      
      for (r = 0; r < length; r++) {
         fprintf(file_val, "%-3d  %-3d  %-3d  %-5.1lf  %-5.1lf  %-3d  %+.15lf  %+.15lf  %+.15lf  %+.15lf\n",
                 site_start, site_start + r, r,
                 (2.0*(double)site_start - 1.0)*0.5, (2.0*site_start + 2.0*r + 1.0)*0.5, r + 1,
                 Val_Real_Right[r][i], Val_Img_Right[r][i], Val_Real_Left[r][i], Val_Img_Left[r][i]);
      }
      fprintf(file_val, "\n");
      fclose(file_val);
   }
   
   for (i = 0; i < out_num; i++) {
      sprintf(Vec_D, "./result/%s/%s/Vector_%d", D1, D2, i);
      mkdir(Vec_D, 0777);
      for (r = 0; r < length; r++) {
         sprintf(File_Name_Vec, "%s/%s_Vec_%d_%03d.txt", Vec_D, File_Name, i, r);
         file_vec = fopen(File_Name_Vec, "a+");
         fprintf(file_vec,"###N=%d,Ne=%d,LocSpin=%1.1lf,P=%d,BC=%s\n###t=%.1lf,J=%.3lf,I_xy=%.3lf,I_z=%.3lf,D_z=%.3lf,h_x=%.3lf,mu=%.5lf\n",
                 Model->tot_site,
                 Model->tot_ele,
                 (double)Model->spin_loc/2.0,
                 Model->tot_parity,
                 Model->BC,
                 Model->t,
                 Model->J,
                 Model->I_xy,
                 Model->I_z,
                 Model->D_z,
                 Model->h_x,
                 Model->mu
                 );
         
         for (j = 0; j < dim; j++) {
            fprintf(file_vec, "%-3d  %-3d  %-3d  %-5.1lf  %-5.1lf  %-3d  %-3d  %+.15lf  %+.15lf  %+.15lf  %+.15lf\n",
                    site_start, site_start + r, r,
                    (2.0*(double)site_start - 1.0)*0.5, (2.0*site_start + 2.0*r + 1.0)*0.5, r + 1,
                    j, Vec_Real_Right[r][i][j], Vec_Img_Right[r][i][j], Vec_Real_Left[r][i][j], Vec_Img_Left[r][i][j]);
         }
         fprintf(file_vec, "\n");
         fclose(file_vec);
      }
   }
   
   FREE_ARRAY_DOUBLE2(Val_Real_Right, length);
   FREE_ARRAY_DOUBLE2(Val_Img_Right , length);
   FREE_ARRAY_DOUBLE2(Val_Real_Left , length);
   FREE_ARRAY_DOUBLE2(Val_Img_Left  , length);
   FREE_ARRAY_DOUBLE3(Vec_Real_Right, length, dim);
   FREE_ARRAY_DOUBLE3(Vec_Img_Right , length, dim);
   FREE_ARRAY_DOUBLE3(Vec_Real_Left , length, dim);
   FREE_ARRAY_DOUBLE3(Vec_Img_Left  , length, dim);
   
}
