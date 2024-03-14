/* compile with this: mex CXXFLAGS="\$CXXFLAGS -std=c++11" mexExtractAnchorsVars.cpp extractconstraints.cpp */
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
    double **bounds, *anchorgrp, *varsgrp, *tmp_list;     //2d array (no. of bounds * 3) (d_i,d_j,d_ij)    
    
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
    size_t len_anchorgrp = mxGetN(prhs[1]);
    //anchorgrp = (double *)mxCalloc(len_q, sizeof(double));
    anchorgrp = mxGetPr(prhs[1]);
    
    // i/p 3
    size_t len_varsgrp = mxGetN(prhs[2]);
    //anchorgrp = (double *)mxCalloc(len_q, sizeof(double));
    varsgrp = mxGetPr(prhs[2]);    
    
    //-------------get vector<int> d_i, d_j, d_ij ----------------//
    vector<int> d_i, d_j, anchor_atoms, vars_atoms;
    vector<double> d_ij;
    
    for(int x = 0; x < row_bounds; x++){
            d_i.push_back(bounds[x][0]);
            d_j.push_back(bounds[x][1]);
            d_ij.push_back(bounds[x][2]);
    }
    
    anchor_atoms.insert(anchor_atoms.begin(), anchorgrp, anchorgrp+len_anchorgrp);
    vars_atoms.insert(vars_atoms.begin(), varsgrp, varsgrp+len_varsgrp);
    
    // -------------- get the extracted (anchor,vars, dist) cons -----------------------//
    vector<pair<pair<int,int>, double>> finalresult = getExtractedConsAnchorAnchor(d_i, d_j, d_ij, anchor_atoms);
    
    //----------------associate the o/p to plhs[0]------------------//
    plhs[0] = mxCreateDoubleMatrix(finalresult.size(), 3, mxREAL);
    double *anchorsAnchor = mxGetPr(plhs[0]);
    
    //----------------assign finalresult back to matlab-------------//    
    int i=0;
    for(vector<pair<pair<int,int>,double>>::iterator it=finalresult.begin(); it!=finalresult.end(); ++it,++i){
            pair<int,int> a = it->first;
            double b = it->second;
             
          //  mexPrintf("\n(i,j,d_ij): (%d,%d,%f) ",a.first,a.second, b);
            anchorsAnchor[i] = a.first;
            anchorsAnchor[i+ finalresult.size()] = a.second;
            anchorsAnchor[i+ 2*finalresult.size()] = b;            
     }

    // -------------- get the extracted (vars, vars) cons -----------------------//
    finalresult = getExtractedConsVarsVars(d_i, d_j, d_ij, vars_atoms);
    
    //----------------associate the o/p to plhs[1]------------------//
    plhs[1] = mxCreateDoubleMatrix(finalresult.size(), 3, mxREAL);
    double *varsVars = mxGetPr(plhs[1]);
    
    //----------------assign finalresult back to matlab-------------//    
    i=0;
    for(vector<pair<pair<int,int>,double>>::iterator it=finalresult.begin(); it!=finalresult.end(); ++it,++i){
            pair<int,int> a = it->first;
            double b = it->second;
             
          //  mexPrintf("\n(i,j,d_ij): (%d,%d,%f) ",a.first,a.second, b);
            varsVars[i] = a.first;
            varsVars[i+ finalresult.size()] = a.second;
            varsVars[i+ 2*finalresult.size()] = b;            
     }
    
}


