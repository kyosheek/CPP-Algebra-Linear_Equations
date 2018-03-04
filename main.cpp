#include "matrix.h"
#include "utils.h"

using namespace std;

void task1();
void task2();
void task3();

int main() {
    cout.setf(ios::fixed);
    cout.precision(5);
    int c;
    cout << "choose a task: ";
    cin >> c;
    switch(c) {
        case 1:
            task1();
            break;
        case 2:
            task2();
            break;
        case 3:
            task3();
            break;
        default:
            cout << "try again";
            break;
    }
    return 0;
}

void task1() {
    ifstream file;
    file.open("/Users/kyoshee/ClionProjects/kp3_matrix/1.txt");

    double k;
    Matrix a,b,c,d;
    file >> k >> a >> b >> c;
    file.close();

    try {
        b.transpose();
        c.reverse();
        b *= k;
        d = b * c;
        a += d;
        cout << a.det() << endl << a;
    } catch(const char* c) {
        printf("%s\n", c);
    }
}

void task2() {
    ifstream file;
    file.open("/Users/kyoshee/ClionProjects/kp3_matrix/2.txt");
    Matrix a;
    file >> a;
    file.close();
    try {
        vector<vector<double>> x = getBasis(a);
    } catch(const char* c) {
        printf("%s\n", c);
    }
}

void task3() {
    ifstream file;
    file.open("/Users/kyoshee/ClionProjects/kp3_matrix/3.txt");
    Matrix a;
    vector<double> b;
    file >> a;
    b.resize(a.columns);
    for (int i = 0; i < a.columns; ++i) {
        file >> b[i];
    }
    file.close();
    try {
        vector<vector<double>> bas = getBasisNH(a, b);
        vector<double> x = solveLnH(a, b);
        ofstream ofile;
        ofile.setf(ios::fixed);
        ofile.open("/Users/kyoshee/ClionProjects/kp3_matrix/out.txt");
        ofile << a;
        ofile.close();
        cout << x.size() << " " << bas.size() << endl;
        for (auto i : x) {
            cout << i << " ";
        }
        cout << endl;
        for (auto i : bas) {
            for (auto j : i) {
                cout << j << " ";
            }
            cout << endl;
        }
        cout << a;
    } catch(const char* c) {
        printf("%s\n", c);
    }
}