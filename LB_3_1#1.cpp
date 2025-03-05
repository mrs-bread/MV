#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
double D(const vector<vector<double>>& matrix, int n) {
    if (n == 1) return matrix[0][0];
    double det = 0;
    vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
    for (int col = 0; col < n; ++col) {
        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j < col)
                    submatrix[i - 1][j] = matrix[i][j];
                else if (j > col)
                    submatrix[i - 1][j - 1] = matrix[i][j];
            }
        }
        det += pow(-1, col) * matrix[0][col] * D(submatrix, n - 1);
    }
    return det;
}

vector<vector<double>> inverseMatrix(const vector<vector<double>>& matrix, int n) {
    double d = D(matrix, n);
    if (fabs(d) < 1e-9) {
        cout << "Матрица вырождена, обратной матрицы не существует." << endl;
        return {};
    }

    vector<vector<double>> adjugate(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
            for (int row = 0; row < n; ++row) {
                for (int col = 0; col < n; ++col) {
                    if (row != i && col != j) {
                        int newRow = row < i ? row : row - 1;
                        int newCol = col < j ? col : col - 1;
                        submatrix[newRow][newCol] = matrix[row][col];
                    }
                }
            }
            adjugate[j][i] = pow(-1, i + j) * D(submatrix, n - 1);
        }
    }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            adjugate[i][j] /= d;

    return adjugate;
}

double Norma(const vector<vector<double>>& matrix, int n, int m) {
    double norma = 0;
    for (int j = 0; j < m; ++j) {
        double col_sum = 0;
        for (int i = 0; i < n; ++i)
            col_sum += fabs(matrix[i][j]);
        if (col_sum > norma)
            norma = col_sum;
    }
    return norma;
}

vector<double> G(vector<vector<double>> A, vector<double> b, int n) {
    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k)
            if (fabs(A[k][i]) > fabs(A[max_row][i]))
                max_row = k;

        swap(A[i], A[max_row]);
        swap(b[i], b[max_row]);

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

vector<double> Nev(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b, int n) {
    vector<double> r(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            r[i] += A[i][j] * x[j];
        r[i] = b[i] - r[i];
    }
    return r;
}

int main() {
    setlocale(LC_ALL, "RUS");
    vector<vector<double>> A1 = { {1.80, 0.80, 0.64},
                                 {1.00, 0.90, 0.81},
                                 {1.00, 1.10, 1.21} };
    vector<double> b1 = { erf(0.80), erf(0.90), erf(1.10) };

    vector<double> x1 = G(A1, b1, A1.size());
    cout << "Решение системы (задача 1):" << endl;
    for (int i = 0; i < A1.size(); ++i)
        cout << "x[" << i << "] = " << x1[i] << endl;

    double normA1 = Norma(A1, A1.size(), A1.size());
    cout << "Норма матрицы A1: " << normA1 << endl;

    vector<vector<double>> invA1 = inverseMatrix(A1, A1.size());
    if (!invA1.empty()) {
        double normInvA1 = Norma(invA1, A1.size(), A1.size());
        cout << "Норма обратной матрицы A1: " << normInvA1 << endl;
        double condA1 = normA1 * normInvA1;
        cout << "Обусловленность матрицы A1: " << condA1 << endl;
    }

    vector<double> r1 = Nev(A1, x1, b1, A1.size());
    cout << "Невязка:" << endl;
    for (int i = 0; i < A1.size(); ++i)
        cout << "r[" << i << "] = " << r1[i] << endl;

    double sum_x1 = 0;
    for (int i = 0; i < A1.size(); ++i)
        sum_x1 += x1[i];
    cout << "Сумма x1 + x2 + x3: " << sum_x1 << endl;
    cout << "Erf(1.0) = " << erf(1.0) << endl << "Разность = " << abs(sum_x1 - erf(1.0)) << endl;
    return 0;
}
