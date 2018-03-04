#include "matrix.h"
#include "utils.h"
#include <algorithm>

/*
 * Triangle matrix
 */
// Changes rows so main diagonal's elements are not zeroes
Matrix changeRows(Matrix& Mx) {
    vector<double> tmp;
    vector<int> counts;
    int count = 0, tmpi;
    for (int i = 0; i < Mx.rows; ++i) {
        for (int j = 0; j < Mx.rows; ++j) {
            if (Mx.elements[i][j] > -0.0000001 && Mx.elements[i][j] < 0.0000001) {
                ++count;
            }
            if (Mx.elements[i][j] <= -0.0000001 || Mx.elements[i][j] >= 0.0000001) {
                break;
            }
        }
        counts.push_back(count);
        count = 0;
    }
    for (int i = 0; i < counts.size() - 1; ++i) // i - номер прохода
    {
        for (int j = 0; j < counts.size() - 1; ++j) // внутренний цикл прохода
        {
            if (counts[j + 1] < counts[j])
            {
                tmp = Mx.elements[j + 1];
                Mx.elements[j + 1] = Mx.elements[j];
                Mx.elements[j] = tmp;
                tmpi = counts[j+1];
                counts[j+1] = counts[j];
                counts[j] = tmpi;
            }
        }
    }
    return Mx;
}
Matrix changeRows(Matrix& Mx, vector<double> &b) {
    vector<double> tmp;
    vector<int> counts;
    int count = 0, tmpi;
    double tmpb;
    for (int i = 0; i < Mx.rows; ++i) {
        for (int j = 0; j < Mx.rows; ++j) {
            if (Mx.elements[i][j] > -0.0000001 && Mx.elements[i][j] < 0.0000001) {
                ++count;
            }
            if (Mx.elements[i][j] != 0) {
                break;
            }
        }
        counts.push_back(count);
        count = 0;
    }
    for (int i = 0; i < counts.size() - 1; ++i) // i - номер прохода
    {
        for (int j = 0; j < counts.size() - 1; ++j) // внутренний цикл прохода
        {
            if (counts[j + 1] < counts[j])
            {
                tmp = Mx.elements[j + 1]; Mx.elements[j + 1] = Mx.elements[j]; Mx.elements[j] = tmp;
                tmpi = counts[j+1]; counts[j+1] = counts[j]; counts[j] = tmpi;
                tmpb = b[j+1]; b[j+1] = b[j]; b[j] = tmpb;
            }
        }
    }
    return Mx;
}

// Eliminates elements under main diagonal
Matrix doColumns(Matrix& Mx) {
    int step = 0;
    double mp = 0, mp2 = 0;
    for (auto row = 0; row != Mx.elements.size(); ++row) {
        for (auto col = 0; col != Mx.elements[0].size(); ++col) {
            if (Mx.elements[row][col] != 0) {
                mp = Mx.elements[row][col];
                step = col;
                break;
            }
        }
        for (auto row2 = row+1; row2 != Mx.elements.size(); ++row2) {
            if (Mx.elements[row2][step] != 0) {
                mp2 = Mx.elements[row2][step] * (1 / mp);
                for (auto col2 = 0; col2 != Mx.elements[0].size(); ++col2) {
                     Mx.elements[row2][col2] -= mp2 * Mx.elements[row][col2];
                }
                Mx.elements[row2][step] = 0; // ROUNDED
            }
        }
    }
    return Mx;
}
Matrix doColumns(Matrix& Mx, vector<double> &b) {
    int step = 0;
    double mp = 0, mp2 = 0;
    for (auto row = 0; row != Mx.elements.size(); ++row) {
        for (auto col = 0; col != Mx.elements[0].size(); ++col) {
            if (Mx.elements[row][col] != 0) {
                mp = Mx.elements[row][col];
                step = col;
                break;
            }
        }
        for (auto row2 = row+1; row2 != Mx.elements.size(); ++row2) {
            if (Mx.elements[row2][step] != 0) {
                mp2 = Mx.elements[row2][step] * (1 / mp);
                for (auto col2 = 0; col2 != Mx.elements[0].size(); ++col2) {
                    Mx.elements[row2][col2] -= mp2 * Mx.elements[row][col2];
                }
                b[row2] -= mp2 * b[row];
                Mx.elements[row2][step] = 0; // ROUNDED
            }
        }
    }
    return Mx;
}
Matrix doColumnsOver(Matrix& Mx, vector<double> &b) {
    int step = 0;
    double mp = 0, mp2 = 0;
    for (auto row = Mx.elements.size()-1; row != -1; --row) {
        for (auto col = Mx.elements[0].size()-1; col != -1; --col) {
            if (Mx.elements[row][col] != 0) {
                mp = Mx.elements[row][col];
                step = col;
                break;
            }
        }
        if (mp != 0) {
            for (auto row2 = row - 1; row2 != -1; --row2) {
                if (Mx.elements[row2][step] != 0) {
                    mp2 = Mx.elements[row2][step] * (1 / mp);
                    for (auto col2 = 0; col2 != Mx.elements[0].size(); ++col2) {
                        Mx.elements[row2][col2] -= mp2 * Mx.elements[row][col2];
                    }
                    b[row2] -= mp2 * b[row];
                    Mx.elements[row2][step] = 0; // ROUNDED
                }
            }
        }
    }
    return Mx;
}

// Makes main diagonal's elements equal to one
Matrix toEDiag(Matrix& Mx, vector<double> &b) {
    double mp = 0;
    for (int row = 0; row < Mx.elements.size(); ++row) {
        for (int col = 0; col < Mx.elements[0].size(); ++col) {
            if (Mx.elements[row][col] != 0) {
                mp = Mx.elements[row][col];
                break;
            }
        }
        for (int col = 0; col < Mx.elements[0].size(); ++col) {
            if (Mx.elements[row][col] != 0) {
                Mx.elements[row][col] *= 1/mp;
            }
        }
        b[row] *= 1/mp;
    }
    return Mx;
}

/*
 * Reverse matrix
 */
// Changes rows
void changeRowsMx(Matrix& Mx, Matrix& E) {
    vector<double> tmpm, tmpe;
    vector<int> counts;
    int count = 0;
    for (int i = 0; i < Mx.rows; ++i) {
        for (int j = 0; j < Mx.rows; ++j) {
            if (Mx.elements[i][j] > -0.0000001 && Mx.elements[i][j] < 0.0000001) {
                ++count;
            }
            if (Mx.elements[i][j] < -0.0000001 || Mx.elements[i][j] > 0.0000001) {
                break;
            }
        }
        counts.push_back(count);
        count = 0;
    }
    for (int i = 0; i < counts.size() - 1; ++i) // i - номер прохода
    {
        for (int j = 0; j < counts.size() - 1; ++j) // внутренний цикл прохода
        {
            if (counts[j + 1] < counts[j])
            {
                tmpm = Mx.elements[j + 1];
                Mx.elements[j + 1] = Mx.elements[j];
                Mx.elements[j] = tmpm;
                tmpe = E.elements[j+1];
                E.elements[j+1] = E.elements[j];
                E.elements[j] = tmpe;
            }
        }
    }
}

void doColumnsOverMx(Matrix& Mx, Matrix& E) {
    double mp = 0, mp2 = 0;
    for (int i = Mx.rows - 1; i != -1; --i) {
        mp = Mx.elements[i][i];
        for (int j = i - 1; j != -1; --j) {
            mp2 = Mx.elements[j][i] / mp;
            for (int k = 0; k < Mx.columns; ++k) {
                Mx.elements[j][k] -= mp2 * Mx.elements[i][k];
                E.elements[j][k] -= mp2 * E.elements[i][k];
            }
        }
    }
}

// Eliminates elements under main diagonal
void doColumnsMx(Matrix& Mx, Matrix& E) {
    double mp = 0;
    double mp2 = 0;
    for (int i = 0; i < Mx.rows; ++i) {
        mp = Mx.elements[i][i];
        for (int j = i + 1; j < Mx.rows; ++j) {
            mp2 = Mx.elements[j][i] / mp;
            for (int k = 0; k < Mx.columns; ++k) {
                Mx.elements[j][k] -= mp2 * Mx.elements[i][k];
                E.elements[j][k] -= mp2 * E.elements[i][k];
            }
        }
    }
}

// Makes main diagonal's elements equal to one
void toEDiag(Matrix& Mx, Matrix& E) {
    double mp = 0;
    for (int row = 0; row < Mx.rows; ++row) {
        mp = Mx.elements[row][row];
        for (int j = 0; j < Mx.columns; ++j) {
            Mx.elements[row][j] /= mp;
            E.elements[row][j] /= mp;
        }
    }
}

// Eliminates elements above main diagonal
void toEMx(Matrix& Mx, Matrix& E) {
    int step = Mx.columns - 1;
    double mp;
    for (int row = Mx.rows - 1; row != -1; --row) {
        if (Mx.elements[row][step] >= 0.9 && Mx.elements[row][step] <= 1.000001) {
            for (int row2 = row - 1; row2 != -1; --row2) {
                if (Mx.elements[row2][step] != 0) {
                    mp = Mx.elements[row2][step];
                    for (int col2 = Mx.columns - 1; col2 != -1; --col2) {
                        Mx.elements[row2][col2] *= 1/mp;
                        E.elements[row2][col2] *= 1/mp;
                    }
                }
            }
        }
        --step;
    }
}

// Solve system of linear homogeneous equations after Gaussian method
Matrix doColumnsE(Matrix& Mx) {
    double mp = 0, mp2 = 0;
    int pos = 0;
    for (int i = 0; i < Mx.rows; ++i) {
        for (int j = 0; j < Mx.columns; ++j) {
            if (Mx.elements[i][j] != 0) {
                mp = Mx.elements[i][j];
                break;
            }
        }
        for (int k = 0; k < Mx.columns; ++k) {
            Mx.elements[i][k] /= mp;
        }
    }
    mp = 0;
    for (int i = Mx.rows - 1; i != -1; --i) {
        for (int j = 0; j < Mx.columns; ++j) {
            if (Mx.elements[i][j] >= 0.9999998 && Mx.elements[i][j] <= 1.0000001) {
                mp = Mx.elements[i][j];
                pos = j;
                break;
            }
        }
        if (mp != 0) {
            for (int k = i - 1; k != -1; --k) {
                    if (Mx.elements[k][pos] <= -0.0000001 || Mx.elements[k][pos] >= 0.0000001) {
                        mp2 = Mx.elements[k][pos] / mp;
                        for (int h = 0; h < Mx.columns; ++h) {
                            Mx.elements[k][h] -= mp2 * Mx.elements[i][h];
                        }
                }
            }
        }
    }

    return Mx;
}
vector<int> getFreeColumns(Matrix& Mx) {
    int count = 0, count2 = 0;
    bool ise = false;
    vector<int> pp;
    for (int i = 0; i < Mx.columns; ++i) {
        for (int j = 0; j < Mx.rows; ++j) {
            if (Mx.elements[j][i] == 1 && find(pp.begin(), pp.end(), j) == pp.end()) ise = true;
            if (Mx.elements[j][i] != 0) ++count;
        }
        if (count == 1 && ise) pp.push_back(i);
        count = 0;
        ise = false;
        ++count2;
    }
    return pp;
}
vector<vector<double>> getBasis(Matrix& a) {
    a.toMainLedge();
    doColumnsE(a);
    //vector<int> freec = getFreeColumns(a);
    vector<vector<double>> x;
    vector<int> freedom;
    int count = 0, size = 0;
    for (int i = 0; i < a.rows; ++i) {
        ++count;
        if (all_of(a.elements[i].begin(), a.elements[i].end(), [](double c){return c == 0;})) {
            ++size;
            freedom.push_back(i);
        }
    }
    if (count < a.columns) {
        size += a.columns - count;
        for (int i = count; i < a.columns; ++i) {
            freedom.push_back(i);
        }
    }

    x.resize(size);
    for (int i = 0; i < x.size(); ++i) {
        x[i].resize(a.columns);
    }
    int step = 0;
    for (int i = 0; i < a.columns; ++i) {
        if (find(freedom.begin(), freedom.end(), i) != freedom.end()) {
            x[step][i] = 1;
            ++step;
        }
    }
    for (int xx = 0; xx < x.size(); ++xx) {
        for (int i = 0; i < a.rows; ++i) {
            for (int j = 0; j < a.columns; ++j) {
                x[xx][i] -= a.elements[i][j] * x[xx][j];
            }
        }
    }
    cout << "!!!basis!!!" << endl;
    cout << x.size() << " " << a.columns << endl;
    for (auto i : x) {
        for (auto j : i) {
            cout << j << " ";
        }
        cout << endl;
    }
    cout << "!!!basis!!!"<< endl;
    return x;
}
void doColumnsE(Matrix& Mx, vector<double>&b) {
    double mp = 0, mp2 = 0;
    int pos = 0;
    for (int i = 0; i < Mx.rows; ++i) {
        for (int j = 0; j < Mx.columns; ++j) {
            if (Mx.elements[i][j] != 0) {
                mp = Mx.elements[i][j];
                break;
            }
        }
        for (int k = 0; k < Mx.columns; ++k) {
            Mx.elements[i][k] /= mp;
        }
        b[i] /= mp;
    }
    mp = 0;
    for (int i = Mx.rows - 1; i != -1; --i) {
        for (int j = 0; j < Mx.columns; ++j) {
            if (Mx.elements[i][j] >= 0.9999998 && Mx.elements[i][j] <= 1.0000001) {
                mp = Mx.elements[i][j];
                pos = j;
                break;
            }
        }
        if (mp != 0) {
            for (int k = i - 1; k != -1; --k) {
                if (Mx.elements[k][pos] <= -0.0000001 || Mx.elements[k][pos] >= 0.0000001) {
                    mp2 = Mx.elements[k][pos] / mp;
                    for (int h = 0; h < Mx.columns; ++h) {
                        Mx.elements[k][h] -= mp2 * Mx.elements[i][h];
                    }
                    b[k] -= mp2 * b[i];
                }
            }
        }
    }
}

vector<vector<double>> getBasisNH(Matrix& a, vector<double>& b) {
    a.toMainLedge(b);
    doColumnsE(a, b);
    //vector<int> freec = getFreeColumns(a);
    vector<vector<double>> x;
    vector<int> freedom;
    int count = 0, size = 0, expandedsize = 0;
    vector<vector<double>> tocopy;
    tocopy.resize(a.rows);
    for (int i = 0; i < a.rows; ++i) {
        tocopy[i].resize(a.columns+1);
        for (int j = 0; j < a.columns; ++j) {
            tocopy[i][j] = a.elements[i][j];
        }
        tocopy[i][a.columns] = b[i];
    }
    Matrix expanded(a.rows, a.columns+1, tocopy);
    for (int i = 0; i < a.rows; ++i) {
        ++count;
        if (all_of(a.elements[i].begin(), a.elements[i].end(), [](double c){return c > -0.0000001 && c < 0.0000001;})) {
            ++size;
            freedom.push_back(i);
        }
    }
    for (int i = 0; i < expanded.rows; ++i) {
        ++count;
        if (all_of(expanded.elements[i].begin(), expanded.elements[i].end(), [](double c){return c > -0.0000001 && c < 0.0000001;})) {
            ++expandedsize;
        }
    }
    /*
    * teorema Kronkera - Kapelli
    */
    if (expandedsize != size) {
        throw "rank of matrix is not equal to rank of expanded matrix";
    }
    if (count < a.columns) {
        size += a.columns - count;
        for (int i = count; i < a.columns; ++i) {
            freedom.push_back(i);
        }
    }

    x.resize(size);
    for (int i = 0; i < x.size(); ++i) {
        x[i].resize(a.columns);
    }
    int step = 0;
    for (int i = 0; i < a.columns; ++i) {
        if (find(freedom.begin(), freedom.end(), i) != freedom.end()) {
            x[step][i] = 1;
            ++step;
        }
    }
    for (int xx = 0; xx < x.size(); ++xx) {
        for (int i = 0; i < a.rows; ++i) {
            for (int j = 0; j < a.columns; ++j) {
                x[xx][i] -= a.elements[i][j] * x[xx][j];
            }
        }
    }
    cout << "!!!basis!!!"<< endl;
    cout << x.size() << " " << a.columns << endl;
    for (auto i : x) {
        for (auto j : i) {
            cout << j << " ";
        }
        cout << endl;
    }
    cout << "!!!basis!!!"<< endl;
    return x;
}

// Solve system of linear non-homogeneous equations after Gaussian method
vector<double> solveLnH(Matrix& a, vector<double> &b) {
    vector<double> x;
    x.resize(a.columns);
    for (int i = a.rows - 1; i > -1; --i) {
        if (all_of(a.elements[i].begin(), a.elements[i].end(), [](int c){return c == 0;})) {
            x[i] = 1;
        } else {
            x[i] = b[i];
            for (int j = 0; j < a.columns; ++j) {
                if (a.at(i,j) != 0) {
                    for (int k = j+1; k < a.columns; ++k) {
                        x[i] -= a.at(i, k) * x[k];
                    }
                }
            }
        }
    }
    return x;
}