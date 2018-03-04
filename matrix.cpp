#include "matrix.h"
#include "utils.h"

using namespace std;

/*
 * Constructors
 */
Matrix::Matrix(int row, int col) {
    this->rows = row;
    this->columns = col;
    if (row == col) {
        this->isSquare = true;
    }
    this->elements.resize(row);
    for (int i = 0; i < row; ++i) {
        this->elements[i].resize(col);
    }
    this->change = 0;
}
Matrix::Matrix(int row, int col, double *els) : Matrix(row,col) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            this->elements[i][j] = els[i*row+j];
        }
    }
}
Matrix::Matrix(int row, int col, double **els) : Matrix(row,col) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            this->elements[i][j] = els[i][j];
        }
    }
}
Matrix::Matrix(int row, int col, vector<double> &els) : Matrix(row, col) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            this->elements[i][j] = els[i*row+j];
        }
    }
}
Matrix::Matrix(int row, int col, vector<vector<double>> &els) : Matrix(row, col) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            this->elements[i][j] = els[i][j];
        }
    }
}

/*
 * Operators
 */
ostream& operator<<(ostream& os, const Matrix& Mx) {
    os << Mx.rows << ' ' << Mx.columns << endl;
    for (auto row = Mx.elements.begin(); row != Mx.elements.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            os.ios::width(7);
            os << *col << ' ';
        }
        os << endl;
    }
    return os;
}
istream& operator>>(istream& is, Matrix& Mx) {
    vector<double> push;
    is >> Mx.rows >> Mx.columns;
    push.resize(Mx.columns);
    for (int i = 0; i < Mx.rows; ++i) {
        for (int j = 0; j < Mx.columns; ++j) {
            is >> push[j];
        }
        Mx.elements.push_back(push);
    }
    if (Mx.rows == Mx.columns) {
        Mx.isSquare = true;
    }
    return is;
}
ifstream& operator>>(ifstream &is, Matrix &Mx) {
    vector<double> push;
    is >> Mx.rows >> Mx.columns;
    push.resize(Mx.columns);
    for (int i = 0; i < Mx.rows; ++i) {
        for (int j = 0; j < Mx.columns; ++j) {
            is >> push[j];
        }
        Mx.elements.push_back(push);
    }
    if (Mx.rows == Mx.columns) {
        Mx.isSquare = true;
    }
    return is;
}
Matrix& Matrix::operator=(const Matrix &Mx) {
    this->rows = Mx.rows;
    this->columns = Mx.columns;
    this->elements = Mx.elements;
    return *this;
}
Matrix& Matrix::operator+=(const Matrix &Mx) {
    try {
        if (this->rows != Mx.rows || this->columns != Mx.columns) {
            throw "different size of matrixes, can accumulate";
        }
    } catch (const char* c) {
        printf("%s\n", c);
    }
    for (int i = 0; i < this->rows; ++i) {
        for (int j = 0; j < this->columns; ++j) {
            this->elements[i][j] += Mx.elements[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator*=(const Matrix &Mx) {
    if (this->rows != Mx.columns) {
        throw "wrong number of rows and columns, cant multiply";
    }
    int j = 0, s = 0, m = 0;
    Matrix newMx(this->rows, Mx.columns);
    for (auto row = newMx.elements.begin(); row != newMx.elements.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            *col = 0;
            for (int i = 0; i < this->columns; ++i) {
                *col += this->elements[j+m][i] * Mx.elements[i][j+s];
            }
            ++s;
        }
        ++m;
        s = 0;
    }
    *this = newMx;
    return *this;
}
Matrix& Matrix::operator*=(double a) {
    for (auto row = this->elements.begin(); row != this->elements.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            *col *= a;
        }
    }
    return *this;
}
Matrix Matrix::operator+(const Matrix& Mx) {
    if (this->rows != Mx.rows || this->columns != Mx.columns) {
        throw "different size of matrixes, can accumulate";
    }
    int i = 0, j = 0;
    Matrix newMx(this->rows, this->columns);
    for (auto row = newMx.elements.begin(); row != newMx.elements.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            newMx.elements[i][j] = *col + Mx.elements[i][j];
            ++j;
        }
        ++i;
    }
    return newMx;
}
Matrix Matrix::operator*(const Matrix& Mx) {
    if (this->columns != Mx.rows) {
        throw "wrong number of rows and columns, cant multiply";
    }
    int j = 0, s = 0, m = 0;
    Matrix newMx(this->rows, Mx.columns);
    for (auto row = newMx.elements.begin(); row != newMx.elements.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            *col = 0;
            for (int i = 0; i < this->columns; ++i) {
                *col += this->elements[j+m][i] * Mx.elements[i][j+s];
            }
            ++s;
        }
        ++m;
        s = 0;
    }
    return newMx;
}

/*
 * Methods
 */
double& Matrix::at(int i, int j) {
    if (i > this->rows || j > this->columns) {
        throw "wrong position of element, cant return element";
    }
    return elements[i][j];
}
double Matrix::at(int i, int j) const {
    if (i > this->rows || j > this->columns) {
        throw "wrong position of element, cant return element";
    }
    return elements[i][j];
}
Matrix& Matrix::transpose() {
    Matrix Mx(this->columns, this->rows);
    int i = this->rows-1, j = this->columns-1;
    for (auto row = this->elements.rbegin(); row != this->elements.rend(); ++row) {
        for (auto col = row->rbegin(); col != row->rend(); ++col) {
            Mx.elements[j][i] = *col;
            --j;
        }
        j = this->columns-1;
        --i;
    }
    *this = Mx;
    return *this;
}
Matrix& Matrix::toLedge() {
    changeRows(*this);
    doColumns(*this);
    return *this;
}
Matrix& Matrix::toMainLedge() {
    doColumns(*this);
    changeRows(*this);
    //toEDiag(*this);
    return *this;
}
Matrix& Matrix::toMainLedge(vector<double> &b) {
    doColumns(*this, b);
    changeRows(*this, b);
    toEDiag(*this, b);
    return *this;
}
double Matrix::det() {
    if (!this->isSquare) {
        throw "matrix is not square, cant calculate determinant";
    }
    double det = 1;
    int mp = 0;
    Matrix copy = *this;
    copy.toLedge();
    for (int i = 0; i < copy.rows; ++i) {
        det *= copy.elements[i][i];
    }
    if (copy.change % 2 == 0 || copy.change == 0) { mp = 1; }
    else { mp = -1; }
    det *= mp;
    return det;
}
Matrix& Matrix::reverse() {
    if (this->det() == 0) {
        throw "determinant equals to 0, cant find reverse matrix";
    }
    if (!this->isSquare) {
        throw "matrix is not square, cant find reverse matrix";
    }
    Matrix E(this->rows, this->columns);
    for (int i = 0; i < E.rows; ++i) {
        E.elements[i][i] = 1;
    }
    changeRowsMx(*this, E);
    doColumnsMx(*this, E);
    toEDiag(*this, E);
    doColumnsOverMx(*this, E);
    *this = E;
    return *this;
}

/*
 * Destructor
 */
Matrix::~Matrix() {
    for (auto row : this->elements) {
        row.clear();
    }
    this->elements.clear();
}