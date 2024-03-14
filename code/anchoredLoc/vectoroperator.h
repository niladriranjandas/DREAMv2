// http://www.cplusplus.com/forum/beginner/127918/

#ifndef VECTOROPERATOR_H
#define VECTOROPERATOR_H

#include <vector>
#include <iostream>

using namespace std;

//-----vectors-----------------------------------
vector<double> makevector(double *a, int dim); 

vector<double> operator+(const vector<double>& lhs, const vector<double>& rhs);  //vector +
 
vector<double> operator-(const vector<double>& lhs, const vector<double>& rhs);  //vector -

vector<double> operator*(const double lhs, const vector<double>& rhs);           // scalar multiply vector

double dot(const vector<double>& lhs, const vector<double>& rhs);

double l2norm(vector<double>& arg);

//-----matrices---------------------------------
vector<vector<double> > makeMatrix(double **a, int row, int col);  

vector<vector<double> > operator+(const vector<vector<double> >& lhs, const vector<vector<double> >& rhs); // matrix add

vector<vector<double> > operator+(const vector<vector<double> >& lhs, const vector<vector<double> >& rhs); // matrix sub

vector<vector<double> > operator->*(const double lhs, const vector<vector<double> >& rhs); // scalar multiply with matrix 

#endif 
