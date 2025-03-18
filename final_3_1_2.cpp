#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
using namespace std;
pair<vector<double>, vector<double>> G(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<double> x(n), nev(n);
    for (int k = 0; k < n; k++) {
        int max_row = k;
        for (int i = k; i < n; i++) {
            if (abs(A[i][k]) > abs(A[max_row][k])) {
                max_row = i;
            }
        }
        swap(A[k], A[max_row]);
        swap(b[k], b[max_row]);
        if (abs(A[k][k]) < 1e-12) {
            cerr << "Матрица вырождена!" << endl;
            return {};
        }
        for (int i = k + 1; i < n; i++) {
            double F = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= F * A[k][j];
            }
            b[i] -= F * b[k];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        nev[i] = b[i] - sum;
    }
    return { x, nev };
}
void P_S(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    cout << "Система уравнений:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "  ";
        for (int j = 0; j < n; j++) {
            cout << setw(10) << fixed << setprecision(2) << A[i][j] << (j < n - 1 ? "x" + to_string(j + 1) + " + " : "x" + to_string(j + 1));
        }
        cout << " = " << setw(10) << fixed << setprecision(2) << b[i] << endl;
    }
}
void P_R(const vector<double>& x, const vector<double>& residuals) {
    cout << "\nРешение системы:" << endl;
    for (size_t i = 0; i < x.size(); i++) {
        cout << "  x" << i + 1 << " = " << setw(15) << fixed << setprecision(6) << x[i] << endl;
    }
    cout << "\nНевязки:" << endl;
    for (double res : residuals) {
        cout << "  " << setw(15) << scientific << setprecision(6) << res << endl;
    }
    cout << endl;
}
void A_P(const vector<vector<double>>& A, const vector<double>& b, const string& title) {
    cout << "\n========================================" << endl;
    cout << title << endl;
    P_S(A, b);
    auto result = G(A, b);
    if (!result.first.empty()) {
        P_R(result.first, result.second);
    }
}
int main() {
    setlocale(LC_ALL, "RUS");
    vector<vector<double>> A1 = { {pow(10,-4), 1}, {1, 2} };
    vector<double> b1 = { 1, 4 };
    vector<vector<double>> A2 = {
        {2.34, -4.21, -11.61},
        {8.04, 5.22, 0.27},
        {3.92, -7.99, 8.37}
    };
    vector<double> b2 = { 14.41, -6.44, 55.56 };
    vector<vector<double>> A3 = {
        {4.43, -7.21, 8.05, 1.23, -2.56},
        {-1.29, 6.47, 2.96, 3.22, 6.12},
        {6.12, 8.31, 9.41, 1.78, -2.88},
        {-2.57, 6.93, -3.74, 7.41, 5.55},
        {1.46, 3.62, 7.83, 6.25, -2.35}
    };
    vector<double> b3 = { 2.62, -3.97, -9.12, 8.11, 7.23 };
    A_P(A1, b1, "Пример а)");
    A_P(A2, b2, "Пример б)");
    A_P(A3, b3, "Пример в)");
    return 0;
}
