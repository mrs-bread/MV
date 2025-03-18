//1
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
//2
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;
int rankk(vector<vector<double>>& matrix) {
    int rank = 0;
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<int> basis_col(rows, -1);
    int col = 0;
    for (int row = 0; row < rows && col < cols; ++row) {
        int max_row = row;
        for (int i = row + 1; i < rows; ++i) {
            if (abs(matrix[i][col]) > abs(matrix[max_row][col])) {
                max_row = i;
            }
        }
        if (abs(matrix[max_row][col]) < 1e-6) {
            row--;
            col++;
            if (col >= cols) break;
            continue;
        }
        if (max_row != row) {
            swap(matrix[row], matrix[max_row]);
        }
        basis_col[row] = col;
        rank++;
        double basis = matrix[row][col];
        for (int j = col; j < cols; ++j) {
            matrix[row][j] /= basis;
        }
        for (int i = 0; i < rows; ++i) {
            if (i != row) {
                double factor = matrix[i][col];
                for (int j = col; j < cols; ++j) {
                    matrix[i][j] -= factor * matrix[row][j];
                }
            }
        }
        col++;
    }
    return rank;
}

int main() {
    setlocale(LC_ALL, "RUS");
    vector<vector<double>> A = {
        {0.1, 0.2, 0.3},
        {0.4, 0.5, 0.6},
        {0.7, 0.8, 0.9}};
    vector<double> b = { 0.1, 0.3, 0.5 };
    vector<vector<double>> BigM = A;
    for (size_t i = 0; i < b.size(); ++i) {
        BigM[i].push_back(b[i]);
    }
    vector<vector<double>> A_copy = A;
    vector<vector<double>> BigM_copy = BigM;
    int rankA = rankk(A_copy);
    int rankBig = rankk(BigM_copy);
    cout << "Ранг матрицы A: " << rankA << endl;
    cout << "Ранг расширенной матрицы [A|b]: " << rankBig << endl;
    if (rankA == rankBig && rankA < A[0].size()) {
        cout << "Система имеет множество решений." << endl;
        int n = A[0].size();
        int r = rankA;
        vector<int> basis_cols;
        vector<int> free_cols;
        vector<vector<double>> step_matrix = BigM;
        rankk(step_matrix);
        for (int i = 0; i < n; ++i) {
            bool is_basis = false;
            for (int j = 0; j < r; ++j) {
                bool found = false;
                for (int k = 0; k < step_matrix.size(); ++k) {
                    if (abs(step_matrix[k][i] - 1.0) < 1e-6) {
                        found = true;
                        break;
                    }
                }
                if (found) {
                    basis_cols.push_back(i);
                    is_basis = true;
                    break;
                }
            }
            if (!is_basis) {
                free_cols.push_back(i);
            }
        }
        cout << "Базисные переменные: ";
        for (int col : basis_cols) {
            cout << "x" << col + 1 << " ";
        }
        cout << endl << "Свободные переменные: ";
        for (int col : free_cols) {
            cout << "x" << col + 1 << " ";
        }
        cout << endl << "Общее решение:" << endl;
        for (int i = 0; i < r; ++i) {
            int basis_col_index = -1;
            for (int j = 0; j < n; ++j) {
                bool found = false;
                for (int k = 0; k < step_matrix.size(); ++k) {
                    if (abs(step_matrix[k][j] - 1.0) < 1e-6) {
                        found = true;
                        basis_col_index = j;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            cout << "x" << basis_col_index + 1 << " = ";
            bool first_term = true;
            for (int j = 0; j < free_cols.size(); ++j) {
                double coeff = -step_matrix[i][free_cols[j]];
                if (abs(coeff) > 1e-6) {
                    if (!first_term) {
                        if (coeff > 0)
                            cout << " + ";
                        else
                            cout << " - ";
                    }
                    else if (coeff < 0) {
                        cout << "-";
                    }
                    cout << abs(coeff) << "*x" << free_cols[j] + 1;
                    first_term = false;
                }
            }
            double constant = step_matrix[i][n];
            if (abs(constant) > 1e-6) {
                if (!first_term) {
                    if (constant > 0)
                        cout << " + ";
                    else
                        cout << " - ";
                }
                else cout << "";
                cout << fixed << setprecision(2) << abs(constant);
                if (constant < 0)
                    cout << "";
            }
            cout << endl;
        }
    }
    else if (rankA == rankBig) {
        cout << "Система имеет единственное решение." << endl;}
    else {
        cout << "Система не имеет решений." << endl;}
    return 0;
}
