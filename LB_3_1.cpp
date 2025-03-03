
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Функция для вычисления матричной нормы
double matrixNorm(const vector<vector<double>>& A) {
    double norm = 0;
    for (const auto& row : A) {
        double sum = 0;
        for (double val : row) sum += abs(val);
        norm = max(norm, sum);
    }
    return norm;
}

// Функция вычисления обратной матрицы
vector<vector<double>> inverseMatrix(vector<vector<double>> A) {
    int n = A.size();
    vector<vector<double>> inv(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) inv[i][i] = 1;

    for (int i = 0; i < n; i++) {
        double diag = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] /= diag;
            inv[i][j] /= diag;
        }
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; j++) {
                    A[k][j] -= factor * A[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }
    return inv;
}

// Умножение матрицы на вектор
vector<double> multiplyMatrixVector(const vector<vector<double>>& A, const vector<double>& x) {
    int n = A.size();
    vector<double> result(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

// Разница векторов
vector<double> vectorDifference(const vector<double>& a, const vector<double>& b) {
    int n = a.size();
    vector<double> diff(n);
    for (int i = 0; i < n; i++) diff[i] = a[i] - b[i];
    return diff;
}

// Функция для вычисления определителя матрицы
double determinant(vector<vector<double>> A) {
    int n = A.size();
    double det = 1;
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[pivot][i])) pivot = j;
        }
        if (abs(A[pivot][i]) < 1e-9) return 0;
        if (i != pivot) swap(A[i], A[pivot]), det *= -1;
        det *= A[i][i];
        for (int j = i + 1; j < n; j++) {
            A[i][j] /= A[i][i];
        }
        for (int j = i + 1; j < n; j++) {
            for (int k = i + 1; k < n; k++) {
                A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }
    return det;
}

// Функция решения системы уравнений
void solveSystem(vector<vector<double>> A, vector<double> B) {
    double detA = determinant(A);
    if (abs(detA) < 1e-9) {
        cout << "Матрица вырожденная, система имеет множество решений." << endl;
        return;
    }

    vector<vector<double>> A_inv = inverseMatrix(A);
    vector<double> X = multiplyMatrixVector(A_inv, B);

    cout << fixed << setprecision(6);
    cout << "Решение: x1 = " << X[0] << ", x2 = " << X[1] << ", x3 = " << X[2] << endl;

    double cond_A = matrixNorm(A) * matrixNorm(A_inv);
    cout << "Обусловленность матрицы: " << cond_A << endl;

    vector<double> residual = vectorDifference(multiplyMatrixVector(A, X), B);
    cout << "Невязка: (" << residual[0] << ", " << residual[1] << ", " << residual[2] << ")" << endl;

    double sumX = X[0] + X[1] + X[2];
    cout << "Сумма решений: " << sumX << ", erf(1.0): " << erf(1.0) << endl;
}

int main() {
    setlocale(LC_ALL, "RUS");
    vector<vector<double>> A1 = { {0.8, 1.0, 0.8}, {1.0, 0.9, 0.81}, {1.0, 1.1, 1.21} };
    vector<double> B1 = { erf(0.8), erf(1.0), erf(1.1) };

    cout << "Задача 1: " << endl;
    solveSystem(A1, B1);

    cout << "\nЗадача 2: " << endl;
    vector<vector<double>> A2 = { {0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}, {0.7, 0.8, 0.9} };
    vector<double> B2 = { 0.1, 0.3, 0.5 };
    solveSystem(A2, B2);

    return 0;
}
