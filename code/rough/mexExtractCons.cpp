/* compile with this: mex CXXFLAGS="\$CXXFLAGS -std=c++11" mexExtractCons.cpp extractconstraints.cpp */
/* doesn't works in matlab 2017a but but not in 2014a */
#include<iostream>
#include<vector>
#include<map>

#include "extractconstraints.h"

#include <mex.h>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
   //----------------set up inputs-----------------------------------------//
    double **bounds, *atomgrp, *tmp_list;     //2d array (no. of bounds * 3) (d_i,d_j,d_ij)    
    
    // extract the i/p from matlab call: p = prhs[0]
    size_t col_bounds  = mxGetN(prhs[0]);
    size_t row_bounds  = mxGetM(prhs[0]);
    
    bounds = (double **)mxCalloc(row_bounds, sizeof(double *));
    for(int x=0; x < row_bounds; x++)
        bounds[x] = (double *)mxCalloc(col_bounds,sizeof(double));
    
    tmp_list = mxGetPr(prhs[0]);
    for (int col=0; col < col_bounds; col++) {
         for (int row=0; row < row_bounds; row++) {
              bounds[row][col] = tmp_list[col*row_bounds+row];
          }
    }
    
    // i/p 2
    size_t len_atomgrp = mxGetN(prhs[1]);
    //atomgrp = (double *)mxCalloc(len_q, sizeof(double));
    atomgrp = mxGetPr(prhs[1]);
    
    //-------------get vector<int> d_i, d_j, d_ij ----------------//
    vector<int> d_i, d_j, grp_atoms;
    vector<double> d_ij;
    
    for(int x = 0; x < row_bounds; x++){
            d_i.push_back(bounds[x][0]);
            d_j.push_back(bounds[x][1]);
            d_ij.push_back(bounds[x][2]);
    }
    
    grp_atoms.insert(grp_atoms.begin(), atomgrp, atomgrp+len_atomgrp);
    
    // -------------- get the extracted cons -----------------------//
    vector<pair<pair<int,int>, double>> finalresult = getExtractedCons(d_i, d_j, d_ij, grp_atoms);
    
    //----------------associate the o/p to plhs[0]------------------//
    plhs[0] = mxCreateDoubleMatrix(finalresult.size(), 3, mxREAL);
    double *output = mxGetPr(plhs[0]);
    
    //----------------assign finalresult back to matlab-------------//    
    int i=0;
    for(vector<pair<pair<int,int>,double>>::iterator it=finalresult.begin(); it!=finalresult.end(); ++it,++i){
            pair<int,int> a = it->first;
            double b = it->second;
             
          //  mexPrintf("\n(i,j,d_ij): (%d,%d,%f) ",a.first,a.second, b);
            output[i] = a.first;
            output[i+ finalresult.size()] = a.second;
            output[i+ 2*finalresult.size()] = b;            
     }
    
}


