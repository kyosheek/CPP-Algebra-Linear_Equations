#ifndef KP3_MATRIX_UTILS_H
#define KP3_MATRIX_UTILS_H

#include "matrix.h"

// Gaussian method for upper-triangle matrix
Matrix changeRows(Matrix& Mx);
Matrix changeRows(Matrix& Mx, vector<double> &b);
Matrix doColumns(Matrix& Mx);
Matrix doColumns(Matrix& Mx, vector<double> &b);
Matrix doColumnsOver(Matrix& Mx, vector<double> &b);
Matrix toEDiag(Matrix& Mx, vector<double> &b);

// Methods to inverse matrix
void changeRowsMx(Matrix& Mx, Matrix& E);
void doColumnsMx(Matrix& Mx, Matrix& E);
void doColumnsOverMx(Matrix& Mx, Matrix& E);
void toEDiag(Matrix& Mx, Matrix& E);

// Get basis of system of linear homogeneous equations
Matrix doColumnsE(Matrix& Mx);
vector<vector<double>> getBasis(Matrix& a);

// Solve system of non-homogeneous linear equations
void doColumnsE(Matrix& Mx, vector<double>& b);
vector<vector<double>> getBasisNH(Matrix& a, vector<double>& b);
vector<double> solveLnH(Matrix& a, vector<double> &b);

#endif //KP3_MATRIX_UTILS_H
