#include "vectoroperator.h"
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>

using std::vector;

vector<double> makevector(double *a, int dim){
    vector<double> vect;

     for(int i=0; i<dim; i++)
          vect.push_back(*(a+i));
   return vect;
}

vector<double> operator+(const vector<double>& lhs, const vector<double>& rhs){	// return type is a vector of integers


	if(lhs.size() != rhs.size()){	// Vectors must be the same size in order to add them!
		throw std::runtime_error("Can't add two vectors of different sizes!");
	}

	vector<double> result;	// Declaring the resulting vector, result

	for(int i=0; i < lhs.size(); i++){	// adding each element of the result vector
		result.push_back(lhs.at(i) + rhs.at(i));	// by adding each element of the two together
	}

	return result;	// returning the vector "result"

}

vector<double> operator-(const vector<double>& lhs, const vector<double>& rhs){	// return type is a vector of integers


	if(lhs.size() != rhs.size()){	// Vectors must be the same size in order to add them!
		throw std::runtime_error("Can't add two vectors of different sizes!");
	}

	vector<double> result;	// Declaring the resulting vector, result

	for(int i=0; i < lhs.size(); i++){	// adding each element of the result vector
		result.push_back(lhs.at(i) - rhs.at(i));	// by adding each element of the two together
	}

	return result;	// returning the vector "result"

}

double dot(const vector<double>& lhs, const vector<double>& rhs){

	if(lhs.size() != rhs.size()){	// Vectors must be the same size in order to add them!
		throw std::runtime_error("Can't add two vectors of different sizes!");
	}

      double sum=0;	// Declaring the resulting vector, result

	for(int i=0; i < lhs.size(); i++){	// adding each element of the result vector
	 	sum+= lhs.at(i) * rhs.at(i);	// by adding each element of the two together
	}

	return sum;	// returning the vector "result"

}

vector<double> operator*(const double lhs, const vector<double>& rhs){
  
        vector<double> result;	// Declaring the resulting vector, result

	for(int i=0; i < rhs.size(); i++){	// adding each element of the result vector
		result.push_back(lhs * rhs.at(i));	// by adding each element of the two together
	}

	return result;	// returning the vector "result"
}

double l2norm(vector<double>& arg){
   double ans = 0;
     for(int i=0; i<arg.size(); i++)
         ans+= arg.at(i) * arg.at(i);
   return(sqrt(ans));   
}

//--------------matrix operations---------------------
vector<vector<double> > makeMatrix(double **a, int row, int col){
   vector<vector<double> > result;
   for(int i=0; i < row; i++){
      vector<double> tmp;
      for(int j=0; j < col; j++)
        tmp.push_back(a[i][j]);
      result.push_back(tmp);
   }
   return result;
}

vector<vector<double> > operator+(const vector<vector<double> >& lhs, const vector<vector<double> >& rhs){
    vector<vector<double> > matrix;
    vector<double> matrix1_row;
    vector<double> matrix2_row;

     for(int i=0; i < lhs.size(); i++){
        matrix1_row = lhs.at(i);
        matrix2_row = rhs.at(i);

        vector<double> res_row;
        for(int j=0; j < matrix1_row.size(); j++)
              res_row.push_back(matrix1_row.at(j) + matrix2_row.at(j));
        matrix.push_back(res_row);
     }
    return matrix;
}

vector<vector<double> > operator-(const vector<vector<double> >& lhs, const vector<vector<double> >& rhs){
    vector<vector<double> > matrix;
    vector<double> matrix1_row;
    vector<double> matrix2_row;

     for(int i=0; i < lhs.size(); i++){
        matrix1_row = lhs.at(i);
        matrix2_row = rhs.at(i);

        vector<double> res_row;
        for(int j=0; j < matrix1_row.size(); j++)
              res_row.push_back(matrix1_row.at(j) - matrix2_row.at(j));
        matrix.push_back(res_row);
     }
    return matrix;
}

vector<vector<double> > operator->*(const double lhs, const vector<vector<double> >& rhs){
    vector<vector<double> > matrix;
    vector<double> matrix2_row;

     for(int i=0; i < rhs.size(); i++){
        matrix2_row = rhs.at(i);

        vector<double> res_row;
        for(int j=0; j < matrix2_row.size(); j++)
              res_row.push_back(lhs * matrix2_row.at(j));
        matrix.push_back(res_row);
     }
    return matrix;
}
