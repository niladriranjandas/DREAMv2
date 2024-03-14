// compile by g++ -std=c++11 testExtractConstraints.cpp extractconstraints.cpp //

#include<iostream>
#include<vector>
#include<map>

#include "extractconstraints.h"

using namespace std;

int main(){
    
    vector<int> dist_i, dist_j, grp_atoms;
    vector<double> dist_ij;
   
    dist_i = {1,1,2,2,3,3,4,4,4,5,6,7,10,10};
    dist_j = {2,3,4,5,5,6,6,7,8,1,2,3,5,8};
    dist_ij = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    grp_atoms = {1,3,5,10};

    vector<pair<pair<int,int>, double>> finalresult = getExtractedCons(dist_i, dist_j, dist_ij, grp_atoms);

    for(vector<pair<pair<int,int>,double>>::iterator it=finalresult.begin(); it!=finalresult.end(); ++it){
           pair<int,int> a = it->first;
           double b = it->second;

           printf("\n(i,j,d_ij): (%d,%d,%f) ",a.first,a.second, b);
    }
    cout<<endl;

  //
  vector<int> anchors, vars;
  anchors = {1,4,10};
  vars = {1,3,5,10};

  finalresult = getExtractedConsAnchorAnchor(dist_i, dist_j, dist_ij, anchors);
    for(vector<pair<pair<int,int>,double>>::iterator it=finalresult.begin(); it!=finalresult.end(); ++it){
           pair<int,int> a = it->first;
           double b = it->second;

           printf("\n(i,j,d_ij): (%d,%d,%f) ",a.first,a.second, b);
    }
    cout<<endl;

 finalresult = getExtractedConsVarsVars(dist_i, dist_j, dist_ij, vars);
    for(vector<pair<pair<int,int>,double>>::iterator it=finalresult.begin(); it!=finalresult.end(); ++it){
           pair<int,int> a = it->first;
           double b = it->second;

           printf("\n(i,j,d_ij): (%d,%d,%f) ",a.first,a.second, b);
    }
    cout<<endl;
}
