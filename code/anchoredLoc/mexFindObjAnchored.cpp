/* ******************************************************************************
    Function which returns the objective function (with anchor)
     
    Input:  prhs[0] = X0: initial point (dim * No. of vars)
            prhs[1] = A : anchors (dim * No. anchors)

            prhs[2] = w_ex: (i,j,ex) equality bound between vars
            prhs[3] = w_ea: (i,j,ea) equality bound between var-i and anchor-j

            prhs[4] = w_ux: (i,j,ux) upper bound between vars
            prhs[5] = w_ua: (i,j,ua) upper bound between var-i and anchor-j

            prhs[6] = w_lx: (i,j,lx) lower bound between vars
            prhs[7] = w_la: (i,j,la) lower bound between var-i and anchor-j

            prhs[8] = w[1) : weight for equality bound between vars
                      w[2) : weight for equality bound between var and anchor
                      w[3) : weight for upper bound between vars
                      w[4) : weight for upper bound between var and anchor
                      w[5) : weight for lower bound between vars
                      w[6) : weight for lower bound between var and anchor
                      w[7) : weight for regularization term

   ****************************************************************************
*/

#include "mex.h"
#include "vectoroperator.h"
#include <vector>
#include <cmath>

using std::vector;

double **get2Darray(const mxArray *ptr, int rows, int cols){         
   double **res = (double **)mxCalloc(rows, sizeof(double *));  
   double *tmp_list;
   
   for(int x = 0; x < rows; x++) {
       res[x] = (double *) mxCalloc(cols, sizeof(double));
   }
   tmp_list = mxGetPr(ptr);
   for (int col=0; col < cols; col++) {
       for (int row=0; row < rows; row++) {
           res[row][col] = tmp_list[col*rows+row];
       }
   }
 return res;
}

double *getColi(double **matrix,int row, int col, int pos){
  double *res=(double *)mxCalloc(row,sizeof(double));
      for(int i=0; i < row; i++)
           res[i] = matrix[i][pos];
  return res;
}

void myFreeArray(double **arr,int row_arr,int col_arr){
     for(int i=0;i<row_arr;i++)
         mxFree(arr[i]);
     mxFree(arr);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

   double **X0, *var_indx, **w_ex, **w_ea, **w_lx, **w_la, **w_ux, **w_ua, *w, *ans;

  //extract  i/p from matlab call-------------------------------------------- 
   size_t col_X0    = mxGetN(prhs[0]);
   size_t row_X0    = mxGetM(prhs[0]);

   //size_t col_A    = mxGetN(prhs[1]);
   //size_t row_A    = mxGetM(prhs[1]); 
   size_t n_var       = mxGetN(prhs[1]);

   size_t col_wex  = mxGetN(prhs[2]);
   size_t row_wex  = mxGetM(prhs[2]);

   size_t col_wea  = mxGetN(prhs[3]);
   size_t row_wea  = mxGetM(prhs[3]);

   size_t col_wux  = mxGetN(prhs[4]);
   size_t row_wux  = mxGetM(prhs[4]); 

   size_t col_wua  = mxGetN(prhs[5]);
   size_t row_wua  = mxGetM(prhs[5]); 

   size_t col_wlx  = mxGetN(prhs[6]);
   size_t row_wlx  = mxGetM(prhs[6]);

   size_t col_wla  = mxGetN(prhs[7]);
   size_t row_wla  = mxGetM(prhs[7]); 

   size_t n_w      = mxGetN(prhs[8]);

      //------------------------validation of i/p------------------------------------------------------//
                if(col_X0 < row_X0)
                   mexErrMsgTxt("Dimenssion of vars matrix is dim * no. of vars.");
                
                if(col_wex > 3)
                    mexErrMsgTxt("Dimenssion of var-var equality should be no. of equality bounds * 3.");
                if(col_wea > 3)
                    mexErrMsgTxt("Dimenssion of var-anchor equality should be no. of equality bounds * 3.");

                if(col_wux > 3)
                    mexErrMsgTxt("Dimenssion of var-var upper should be no. of upper bounds * 3.");
                if(col_wua > 3)
                    mexErrMsgTxt("Dimenssion of var-anchor upper should be no. of upper bounds * 3.");

                if(col_wlx > 3)
                    mexErrMsgTxt("Dimenssion of var-var lower should be no. of lower bounds * 3.");
                if(col_wla > 3)
                    mexErrMsgTxt("Dimenssion of var-anchor lower should be no. of lower bounds * 3.");

                if(n_w > 7)
                    mexErrMsgTxt("Dimenssion of weight matrix should be 7 * 1.");

      //-----------------------------------------------------------------------------------------------//
   X0   = get2Darray(prhs[0],row_X0,col_X0);
   //A    = get2Darray(prhs[1],row_A,col_A);
   var_indx    =  mxGetPr(prhs[1]);

   w_ex = get2Darray(prhs[2],row_wex,col_wex);
   w_ea = get2Darray(prhs[3],row_wea,col_wea);

   w_ux = get2Darray(prhs[4],row_wux,col_wux);
   w_ua = get2Darray(prhs[5],row_wua,col_wua);

   w_lx = get2Darray(prhs[6],row_wlx,col_wlx);
   w_la = get2Darray(prhs[7],row_wla,col_wla);
   
   w = mxGetPr(prhs[8]);
   
 // Pointer for the o/p --------------------------------
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    /* pointer to output data */
    ans = mxGetPr(plhs[0]);

 //-----------------------------------------------------
   double tmp, d_ij, dist_xij, dim=3;
   double *col_i;
   vector<double> x_i, x_j, a_j, tmp_;
   *ans = 0;
// mexPrintf("\n grad: X0(row:%f , col:%u)\n",dim,row_X0); 
     // for the equality bound between vars
        tmp = 0;
        for(int i=0; i < row_wex; i++){
           d_ij  = w_ex[i][2];
             col_i = getColi(X0, row_X0, dim, w_ex[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim, w_ex[i][1]-1);
           x_j = makevector(col_i,dim);

           tmp_ = x_i - x_j ;
           tmp += pow(l2norm(tmp_) - d_ij,2);           
        }
        *ans += w[0]*tmp;
        
     // equality bound between var and anchor
        tmp = 0;
        for(int i=0; i < row_wea; i++){
           d_ij  = w_ea[i][2];
             col_i = getColi(X0, row_X0, dim, w_ea[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim, w_ea[i][1]-1);
           a_j = makevector(col_i,dim);
           //a_j = makevector(getColi(A, row_A, dim, w_ea[i][1]),dim);

           tmp_ = x_i - a_j;
           tmp += pow(l2norm(tmp_) - d_ij,2);            
        }
        *ans += w[1] * tmp;

     // upper bound between vars
        tmp = 0;
        for(int i=0; i < row_wux; i++){
           d_ij  = w_ux[i][2];
             col_i = getColi(X0, row_X0, dim, w_ux[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim, w_ux[i][1]-1);
           x_j = makevector(col_i,dim);

           tmp_ = x_i - x_j;
           dist_xij = l2norm(tmp_);

           if(dist_xij > d_ij)
             tmp += pow(dist_xij - d_ij,2);            
        }
        *ans += w[2] * tmp;
        
     // upper bound between vars and anchor
        tmp = 0;
        for(int i=0; i < row_wua; i++){
           d_ij  = w_ua[i][2];
             col_i = getColi(X0, row_X0, dim, w_ua[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim, w_ua[i][1]-1);
           a_j = makevector(col_i,dim);
          //a_j = makevector(getColi(A, row_A, dim, w_ua[i][1]),dim);

           tmp_ = x_i - a_j;
           dist_xij = l2norm(tmp_);

           if(dist_xij > d_ij)
             tmp += pow(dist_xij - d_ij,2);            
           
        }
        *ans += w[3] * tmp;
        
     // lower bound between vars
        tmp = 0;
        for(int i=0; i < row_wlx; i++){
           d_ij  = w_lx[i][2];
             col_i = getColi(X0, row_X0, dim, w_lx[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim, w_lx[i][1]-1);
           x_j = makevector(col_i,dim);

           tmp_ = x_i - x_j;
           dist_xij = l2norm(tmp_);

           if(dist_xij < d_ij) 
               tmp += pow(dist_xij - d_ij,2);            
        }
        *ans += w[4] * tmp;
        
     // lower bound between vars and anchor
        tmp = 0;
        for(int i=0; i < row_wla; i++){
           d_ij  = w_la[i][2];
             col_i = getColi(X0, row_X0, dim, w_la[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim, w_la[i][1]-1);
           a_j = makevector(col_i,dim);
           //a_j = makevector(getColi(A, row_A, dim, w_la[i][1]),dim);

           tmp_ = x_i - a_j;
           dist_xij = l2norm(tmp_);

           if(dist_xij < d_ij) 
             tmp += pow(dist_xij - d_ij,2);            
        }
        *ans += w[5] * tmp;

     // regularization  
        tmp = 0;
        for(int i=0; i < n_var; i++){
             col_i = getColi(X0, row_X0, dim, var_indx[i]-1);
           x_i = makevector(col_i,dim);
           tmp += pow(l2norm(x_i),2);      
        }
        *ans += w[6] * tmp;

        // free memory      
    //--X0---
      myFreeArray(X0,dim,col_X0);
     
    //--w_ex and w_ea--
      myFreeArray(w_ex,row_wex,col_wex);
      myFreeArray(w_ea,row_wea,col_wea);
    // --w_ux and w_ua--
      myFreeArray(w_ux,row_wux,col_wux);
      myFreeArray(w_ua,row_wua,col_wua);
    // -- w_lx and w_la --
      myFreeArray(w_lx,row_wlx,col_wlx);
      myFreeArray(w_la,row_wla,col_wla);
      
    // vector<double> x_i, x_j, a_j, tmp_;
      mxFree(col_i);
      x_i.clear();
      x_j.clear();
      a_j.clear();
      tmp_.clear();
    
       
}

/*

   X0 = (double **)mxCalloc(row_X0, sizeof(double *));  
   for(int x = 0; x < row_X0; x++) {
       X0[x] = (double *) mxCalloc(col_X0, sizeof(double));
   }
   tmp_list = mxGetPr(prhs[0]);
   for (int col=0; col < col_X0; col++) {
       for (int row=0; row < row_X0; row++) {
           X0[row][col] = tmp_list[col*row_X0+row];
      }
   }

   A = (double **)mxCalloc(row_A, sizeof(double *));
   for(int x = 0; x < row_A; x++) {
       A[x] = (double *) mxCalloc(col_A, sizeof(double));
   }
   tmp_list = mxGetPr(prhs[1]);
   for (int col=0; col < col_X0; col++) {
       for (int row=0; row < row_X0; row++) {
           A[row][col] = tmp_list[col*row_X0+row];
      }
   }

   w_ex = (double **)mxCalloc(row_wex, sizeof(double *));
   for(int x = 0; x < row_wex; x++) {
       w_ex[x] = (double *) mxCalloc(col_wex, sizeof(double));
   }   
   tmp_list = mxGetPr(prhs[2]);
   for (int col=0; col < col_wex; col++) {
       for (int row=0; row < row_wex; row++) {
           w_ex[row][col] = tmp_list[col*row_wex+row];
      }
   }

   w_ea = (double **)mxCalloc(row_wea, sizeof(double *));
   for(int x = 0; x < row_wea; x++) {
       w_ea[x] = (double *) mxCalloc(col_wea, sizeof(double));
   }
   tmp_list = mxGetPr(prhs[3]);
   for (int col=0; col < col_wea; col++) {
       for (int row=0; row < row_wea; row++) {
           w_ea[row][col] = tmp_list[col*row_wea+row];
      }
   }

   w_ux = (double **)mxCalloc(row_wux, sizeof(double *));
   for(int x = 0; x < row_wux; x++) {
       w_ux[x] = (double *) mxCalloc(col_wux, sizeof(double));
   }   
   tmp_list = mxGetPr(prhs[4]);
   for (int col=0; col < col_wux; col++) {
       for (int row=0; row < row_wux; row++) {
           w_ux[row][col] = tmp_list[col*row_wux+row];
      }
   }

   w_ua = (double **)mxCalloc(row_wua, sizeof(double *));
   for(int x = 0; x < row_wua; x++) {
       w_ua[x] = (double *) mxCalloc(col_wua, sizeof(double));
   }
   tmp_list = mxGetPr(prhs[5]);
   for (int col=0; col < col_wua; col++) {
       for (int row=0; row < row_wua; row++) {
           w_ua[row][col] = tmp_list[col*row_wux+row];
      }
   }

   w_lx = (double **)mxCalloc(row_wlx, sizeof(double *));
   for(int x = 0; x < row_wlx; x++) {
       w_lx[x] = (double *) mxCalloc(col_wlx, sizeof(double));
   }   
   tmp_list = mxGetPr(prhs[6]);
   for (int col=0; col < col_wlx; col++) {
       for (int row=0; row < row_wlx; row++) {
           w_lx[row][col] = tmp_list[col*row_wlx+row];
      }
   }

   w_la = (double **)mxCalloc(row_wla, sizeof(double *));
   for(int x = 0; x < row_wla; x++) {
       w_la[x] = (double *) mxCalloc(col_wla, sizeof(double));
   }
   tmp_list = mxGetPr(prhs[5]);
   for (int col=0; col < col_wla; col++) {
       for (int row=0; row < row_wla; row++) {
           w_la[row][col] = tmp_list[col*row_wla+row];
      }
   }

*/
