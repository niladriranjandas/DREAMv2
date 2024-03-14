/* ******************************************************************************
    Function which returns the grad of function (with anchor)
     
    Input:  prhs[0] = X0: initial point (dim * No. of vars)
            prhs[1] = A : anchors (dim * No. anchors)

            prhs[2] = w_ex: (i,j,ex) equality bound between vars
            prhs[3] = w_ea: (i,j,ea) equality bound between var-i and anchor-j

            prhs[4] = w_ux: (i,j,ux) upper bound between vars
            prhs[5] = w_ua: (i,j,ua) upper bound between var-i and anchor-j

            prhs[6] = w_lx: (i,j,lx) lower bound between vars
            prhs[7] = w_la: (i,j,la) lower bound between var-i and anchor-j

            prhs[8] = w[0) : weight for equality bound between vars
                      w[1) : weight for equality bound between var and anchor
                      w[2) : weight for upper bound between vars
                      w[3) : weight for upper bound between var and anchor
                      w[4) : weight for lower bound between vars
                      w[5) : weight for lower bound between var and anchor
                      w[6) : weight for regularization term

compile by the following
CFLAGS="$CFLAGS -fopenmp"
CXXFLAGS="$CXXFLAGS -fopenmp"
mex mexFindGradAnchored.cpp vectoroperator.cpp

mex CXXFLAGS="\$CXXFLAGS -std=c++11 -fopenmp" CXXOPTIMFLAGS='\$CXXOPTIMFLAGS -Ofast -DNDEBUG' LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp -O2" mexFindGradAnchored_omp.cpp vectoroperator.cpp
   ****************************************************************************
*/

#include "mex.h"
#include "vectoroperator.h"
#include "omp.h"

using std::vector;

#define NUM_THREADS 6

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

void addAndPlaceAti(double **grad, int row_grad, vector<double> vect, int col_pos){
  if (vect.size() != row_grad) 
     mexErrMsgTxt("incorrect call to addAndPlaceAti: dimenssion mismatch.");

   for(int i=0; i < row_grad; i++)
       grad[i][col_pos] += vect[i];
}

double *getColi(double **matrix,int row, int col, int pos){
  double *res=(double *)mxCalloc(row,sizeof(double));
      for(int i=0; i < row; i++)
           res[i] = matrix[i][pos];
  return res;
}

void placeInColj(double **X, double *a, int row_X, int col_X, int j){
    for(int i=0; i < row_X; i++)
       X[i][j] = *(a+i);
}

void myFreeArray(double **arr,int row_arr,int col_arr){
     for(int i=0;i<row_arr;i++)
         mxFree(arr[i]);
     mxFree(arr);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

   double **X0, *var_indx, **w_ex, **w_ea, **w_lx, **w_la, **w_ux, **w_ua, *w ;

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

      //-------validation of i/p------------------------------------
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
      //------------------------------------------------------------
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
   double dim=3;//dim=col_X0;
    plhs[0] = mxCreateDoubleMatrix(dim, col_X0, mxREAL);
    /* pointer to output data */
    double *ans = mxGetPr(plhs[0]);      
//mexPrintf("\n grad: dim:%f , row_X0:%u, col_X0:%u)\n",dim,row_X0,col_X0);
     int i;
     double **grad = (double **) mxCalloc(dim, sizeof(double *));
        for(i=0; i < dim; i++){
            grad[i] = (double *) mxCalloc(col_X0 ,sizeof(double));
        }    
 //-----------------------------------------------------
   double tmp, d_ij, dist_xij, dist_xiaj, scalar;
   double *col_i;
   vector<double> x_i , x_j, a_j, Xeij, Xiaj;
   
     // for the equality bound between vars
        int n_per_thread;
        omp_set_num_threads(NUM_THREADS);
        omp_set_dynamic(4);
                
        n_per_thread = row_wex/NUM_THREADS;         
        #pragma omp parallel for shared(grad) private(i,d_ij,col_i,x_i,x_j,Xeij,dist_xij,scalar) schedule(static, n_per_thread)
        for(i=0; i < row_wex; i++){
           d_ij  = w_ex[i][2];
             col_i = getColi(X0, row_X0, dim, w_ex[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim, w_ex[i][1]-1);
           x_j = makevector(col_i,dim);

           Xeij = x_i - x_j;
           dist_xij = l2norm(Xeij);
           scalar = 2*w[0]*(1 - d_ij/dist_xij);
          
           addAndPlaceAti(grad, dim, scalar*Xeij, w_ex[i][0]-1);      
           addAndPlaceAti(grad, dim, (-1*scalar)*Xeij, w_ex[i][1]-1);
        }

     // equality bound between var and anchor
        
        n_per_thread = row_wea/NUM_THREADS;
        #pragma omp parallel for shared(grad) private(i,d_ij,col_i,x_i,a_j,Xiaj,dist_xiaj,scalar) num_threads(n_per_thread)
        for(i=0; i < row_wea; i++){
           d_ij  = w_ea[i][2];
             col_i = getColi(X0, row_X0, dim,w_ea[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim,w_ea[i][1]-1);
           a_j = makevector(col_i,dim);
           //a_j = makevector(getColi(A, row_A, dim,w_ea[i][1]),dim);

           Xiaj  = x_i - a_j;
           dist_xiaj = l2norm(Xiaj);
           scalar = 2*w[1]*(1 - d_ij/dist_xiaj);
          
           addAndPlaceAti(grad, dim, scalar*Xiaj, w_ea[i][0]-1);      
        }

// x_i = makevector(getColi(X0, row_X0, dim,w_ux[row_wux-1][0]-1),dim); 
// x_j = makevector(getColi(X0, row_X0, dim,w_ux[row_wux-1][1]-1),dim);
// mexPrintf("\n p: (x_i:%f, x_j:%f, dist:%f",w_ux[row_wux-1][0],w_ux[row_wux-1][1],w_ux[row_wux-1][2]);
// for(int i=0;i<x_i.size();i++)
//     mexPrintf("\t%f",x_i[i]);
// mexPrintf("\n");

     // upper bound between vars
        
        n_per_thread = row_wux/NUM_THREADS;
        #pragma omp parallel for shared(grad) private(i,d_ij,col_i,x_i,x_j,Xeij,dist_xij,scalar) schedule(static, n_per_thread)
        for(i=0; i < row_wux; i++){
           d_ij  = w_ux[i][2];
             col_i = getColi(X0, row_X0, dim,w_ux[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim,w_ux[i][1]-1);
           x_j = makevector(col_i,dim);

           Xeij = x_i - x_j;
           dist_xij = l2norm(Xeij);
           scalar = 2*w[2]*(1 - d_ij/dist_xij);
           
            if(dist_xij > d_ij){
             addAndPlaceAti(grad, dim, scalar*Xeij, w_ux[i][0]-1);      
             addAndPlaceAti(grad, dim, (-1*scalar)*Xeij, w_ux[i][1]-1);
            }
        }

     // upper bound between vars and anchor
        
        n_per_thread = row_wua/NUM_THREADS;
        #pragma omp parallel for shared(grad) private(i,d_ij,col_i,x_i,a_j,Xiaj,dist_xiaj,scalar) schedule(static, n_per_thread)
        for(i=0; i < row_wua; i++){
           d_ij  = w_ua[i][2];
             col_i = getColi(X0, row_X0, dim,w_ua[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim,w_ua[i][1]-1);
           a_j = makevector(col_i,dim);
           //a_j = makevector(getColi(A, row_A, dim,w_ua[i][1]),dim);

           Xiaj = x_i - a_j;
           dist_xiaj = l2norm(Xiaj);
           scalar = 2*w[3]*(1 - d_ij/dist_xiaj);
          
           if (dist_xiaj > d_ij)
                addAndPlaceAti(grad, dim, scalar*Xiaj, w_ua[i][0]-1);
        }

     // lower bound between vars
        
        n_per_thread = row_wlx/NUM_THREADS;
        #pragma omp parallel for shared(grad) private(i,d_ij,col_i,x_i,x_j,Xeij,dist_xij,scalar) schedule(static, n_per_thread)
        for(i=0; i < row_wlx; i++){
           d_ij  = w_lx[i][2];
             col_i = getColi(X0, row_X0, dim,w_lx[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim,w_lx[i][1]-1);
           x_j = makevector(col_i,dim);

           Xeij = x_i - x_j;
           dist_xij = l2norm(Xeij);
           scalar = 2*w[4]*(1 - d_ij/dist_xij);
          
           if(dist_xij < d_ij){
            addAndPlaceAti(grad, dim, scalar*Xeij, w_lx[i][0]-1);      
            addAndPlaceAti(grad, dim, (-1*scalar)*Xeij, w_lx[i][1]-1);
           }
        }

     // lower bound between vars and anchor
        
        n_per_thread = row_wla/NUM_THREADS;
        #pragma omp parallel for shared(grad) private(i,d_ij,col_i,x_i,a_j,Xiaj,dist_xiaj,scalar) schedule(static, n_per_thread)
        for(i=0; i < row_wla; i++){
           d_ij  = w_la[i][2];
             col_i = getColi(X0, row_X0, dim,w_la[i][0]-1);
           x_i = makevector(col_i,dim); 
             col_i = getColi(X0, row_X0, dim,w_la[i][1]-1);
           a_j = makevector(col_i,dim);
           //a_j = makevector(getColi(A, row_A, dim,w_la[i][1]),dim);

           Xiaj = x_i - a_j;
           dist_xiaj = l2norm(Xiaj);
           scalar = 2*w[5]*(1 - d_ij/dist_xiaj);
          
           if (dist_xiaj < d_ij)
                addAndPlaceAti(grad, dim, scalar*Xiaj, w_la[i][0]-1);
        }

     // regularization  
        
        n_per_thread = n_var/NUM_THREADS;
        #pragma omp parallel for shared(grad) private(i,col_i,x_i,scalar) schedule(static, n_per_thread)              
        for(i=0; i < n_var; i++){     
             col_i = getColi(X0, row_X0, dim, var_indx[i]-1);
           x_i = makevector(col_i,dim);
           scalar = 2*w[6];
           addAndPlaceAti(grad, dim, scalar*x_i, var_indx[i]-1); // for the whole matrix X0, they are same, but rowsize reduces l2norm() call
        }

   // assinging "grad" to matlab o/p
        for(int col=0; col < col_X0; col++){
          for(int row=0; row < dim; row++)
             ans[row + col*row_X0] = grad[row][col];
        }
   // free memory      
    //--X0---
      myFreeArray(X0,dim,col_X0);
      
    //--grad---
      myFreeArray(grad,dim,col_X0);   
     
    //--w_ex and w_ea--
      myFreeArray(w_ex,row_wex,col_wex);
      myFreeArray(w_ea,row_wea,col_wea);
    // --w_ux and w_ua--
      myFreeArray(w_ux,row_wux,col_wux);
      myFreeArray(w_ua,row_wua,col_wua);
    // -- w_lx and w_la --
      myFreeArray(w_lx,row_wlx,col_wlx);
      myFreeArray(w_la,row_wla,col_wla);      
      
    //vector<double> x_i , x_j, a_j, Xeij, Xiaj;
      mxFree(col_i);
      x_i.clear();
      x_j.clear();
      a_j.clear();
      Xeij.clear();
      Xiaj.clear();
}

