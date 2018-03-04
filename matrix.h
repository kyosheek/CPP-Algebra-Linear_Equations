#ifndef KP3_MATRIX_MATRIX_H
#define KP3_MATRIX_MATRIX_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

class Matrix{
public:
    int change;
    bool isSquare = false;
    int rows, columns;
    vector<vector<double>> elements;

    /*
     * Constructos
     */
    Matrix() {}
    Matrix(int row, int col);
    Matrix(int row, int col, double *els);
    Matrix(int row, int col, double **els);
    Matrix(int row, int col, vector<double> &els);
    Matrix(int row, int col, vector<vector<double>> &els);
    Matrix(int a) : Matrix(a,a) {};
    Matrix(int a, double *els) : Matrix(a,a,els) {};
    Matrix(int a, double **els) : Matrix(a,a,els) {};
    Matrix(int a, vector<double> &els) : Matrix(a,a,els) {};
    Matrix(int a, vector<vector<double>> &els) : Matrix(a,a,els) {};

    /*
     * Operators
     */
    friend ostream& operator<<(ostream& os, const Matrix &Mx);
    Matrix& operator=(const Matrix& Mx); //copy
    Matrix& operator+=(const Matrix& Mx); //matrix sum assign
    Matrix& operator*=(const Matrix& Mx); //matrix multiply assign
    Matrix& operator*=(double a); //number multiply assign
    Matrix operator+(const Matrix& Mx); //matrix sum
    Matrix operator*(const Matrix& Mx); //matrix multiply

    friend ostream& operator<<(ostream &os, const Matrix &Mx); //output
    friend istream& operator>>(istream &is, Matrix &Mx); //input
    friend ifstream& operator>>(ifstream &is, Matrix &Mx); //input


    /*
     * Methods
     */
    double& at(int i, int j);
    double at(int i, int j) const;
    Matrix& transpose();
    Matrix& toLedge();
    Matrix& toMainLedge();
    Matrix& toMainLedge(vector<double> &b);
    double det();
    Matrix& reverse();

        /*
         * Destructor
         */
    ~Matrix();
};

class sqMatrix : public Matrix {
    sqMatrix() : Matrix(1,1) {}
};

class vMatrix : public Matrix {
    vMatrix() : Matrix(1,1) {}
};

#endif //KP3_MATRIX_MATRIX_H
