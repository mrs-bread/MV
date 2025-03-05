#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
using namespace std;

vector<double> gauss(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    vector<double> x(n, 0.0);//сюда решение

    //это верхний треугольник
    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[max_row][i])) {
                max_row = k;
            }
        }

        if (max_row != i) {
            swap(A[i], A[max_row]);
            swap(b[i], b[max_row]);
        }

        for (int k = i + 1; k < n; ++k) { //строка ниже текующей
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    //решение вычисляем
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) { // справа от диагонального 
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return x;
}

double matrixNorm(const vector<vector<double>>& A) {
    int n = A.size();
    double max_sum = 0.0;
    for (int j = 0; j < n; ++j) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += abs(A[i][j]);
        }
        max_sum = max(max_sum, sum);
    }
    return max_sum;
}


vector<vector<double>> inverseMatrix(vector<vector<double>> A) {
    int n = A.size();
    vector<vector<double>> inv(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        inv[i][i] = 1.0;
    }

    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[max_row][i])) {
                max_row = k;
            }
        }

        if (max_row != i) {
            swap(A[i], A[max_row]);
            swap(inv[i], inv[max_row]);
        }

        double pivot = A[i][i]; //диагональный эл
        for (int j = 0; j < n; ++j) {
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }

    return inv;
}

double obus(const vector<vector<double>>& A) {
    vector<vector<double>> A_copy = A;
    vector<vector<double>> A_inv = inverseMatrix(A_copy);
    return matrixNorm(A) * matrixNorm(A_inv);
}

double nev(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
    int n = A.size();
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        res += abs(sum - b[i]);
    }
    return res;
}

int main() {
    setlocale(LC_ALL, "Russian");
    cout << fixed << setprecision(10);

    vector<vector<double>> A1 = {
        {1.00, 0.80, 0.64},
        {1.00, 0.90, 0.81},
        {1.00, 1.10, 1.21}
    };

    vector<double> b1 = { erf(0.80), erf(0.90), erf(1.10) };

    double cond_A1 = obus(A1);
    cout << "Задача 1:\n";
    cout << "Обусловленность: " << cond_A1 << endl;

    vector<double> x1 = gauss(A1, b1);
    cout << "x1: " << x1[0] << endl;
    cout << "x2: " << x1[1] << endl;
    cout << "x3: " << x1[2] << endl;

    vector<vector<double>> A1_copy = {
        {1.00, 0.80, 0.64},
        {1.00, 0.90, 0.81},
        {1.00, 1.10, 1.21}
    };
    vector<double> b1_copy = { erf(0.80), erf(0.90), erf(1.10) };

    double res1 = nev(A1_copy, x1, b1_copy);
    cout << "Невязка |Ax - b|: " << res1 << endl;

    double sum_x1 = x1[0] + x1[1] + x1[2];
    cout << "x1 + x2 + x3: " << sum_x1 << endl;
    cout << "erf(1.0): " << erf(1.0) << endl;
    // Матрица - это степени 1.1 0,9 и 0,8. Т.е. уравнение параболы вида y = x1 + x2x + x3x^2. Если мы подставим x=1, то получим сумму коэфов. Получается, мы вычисляем значение в точке 1.

    vector<vector<double>> A2 = {
        {0.1, 0.2, 0.3},
        {0.4, 0.5, 0.6},
        {0.7, 0.8, 0.9}
    };

    vector<double> b2 = { 0.1, 0.3, 0.5 };

    cout << "\nЗадача 2:\n";
    cout << "det(A): " << A2[0][0] * (A2[1][1] * A2[2][2] - A2[1][2] * A2[2][1]) - A2[0][1] * (A2[1][0] * A2[2][2] - A2[1][2] * A2[2][0]) + A2[0][2] * (A2[1][0] * A2[2][1] - A2[1][1] * A2[2][0]) << endl;

    return 0;
}
//#include <iostream>
//#include <cmath>
//#include <vector>
//using namespace std;
//
//// Вычисление определителя матрицы
//double determinant(const vector<vector<double>>& A) {
//    int n = A.size();
//    if (n == 1) return A[0][0];
//    if (n == 2) return A[0][0] * A[1][1] - A[0][1] * A[1][0];
//
//    double det = 0;
//    for (int j = 0; j < n; ++j) {
//        vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
//        for (int i = 1; i < n; ++i) {
//            for (int k = 0; k < n; ++k) {
//                if (k < j) submatrix[i - 1][k] = A[i][k];
//                else if (k > j) submatrix[i - 1][k - 1] = A[i][k];
//            }
//        }
//        det += pow(-1, j) * A[0][j] * determinant(submatrix);
//    }
//    return det;
//}
//
//// Вычисление алгебраического дополнения
//double cofactor(const vector<vector<double>>& A, int i, int j) {
//    int n = A.size();
//    vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
//    for (int row = 0, r = 0; row < n; ++row) {
//        if (row == i) continue;
//        for (int col = 0, c = 0; col < n; ++col) {
//            if (col == j) continue;
//            submatrix[r][c++] = A[row][col];
//        }
//        if (row != i) ++r;
//    }
//    return pow(-1, i + j) * determinant(submatrix);
//}
//
//// Вычисление обратной матрицы с помощью алгебраических дополнений
//vector<vector<double>> inverseMatrix(const vector<vector<double>>& A) {
//    int n = A.size();
//    double detA = determinant(A);
//    if (abs(detA) < 1e-9) {
//        cout << "Matrix is singular, cannot compute inverse." << endl;
//        return {};
//    }
//
//    vector<vector<double>> adjugate(n, vector<double>(n));
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            adjugate[j][i] = cofactor(A, i, j);
//        }
//    }
//
//    vector<vector<double>> invA(n, vector<double>(n));
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            invA[i][j] = adjugate[i][j] / detA;
//        }
//    }
//    return invA;
//}
//
//// Реализация метода Гаусса
//void gaussElimination(vector<vector<double>>& A, vector<double>& b) {
//    int n = A.size();
//
//    // Прямой ход метода Гаусса
//    for (int i = 0; i < n; ++i) {
//        // Поиск максимального элемента в столбце
//        int maxRow = i;
//        for (int j = i + 1; j < n; ++j) {
//            if (abs(A[j][i]) > abs(A[maxRow][i])) {
//                maxRow = j;
//            }
//        }
//
//        // Обмен строк
//        swap(A[i], A[maxRow]);
//        swap(b[i], b[maxRow]);
//
//        // Нормализация строки
//        double pivot = A[i][i];
//        for (int j = i; j < n; ++j) {
//            A[i][j] /= pivot;
//        }
//        b[i] /= pivot;
//
//        // Вычитание кратных текущей строки из остальных строк
//        for (int j = i + 1; j < n; ++j) {
//            double factor = A[j][i];
//            for (int k = i; k < n; ++k) {
//                A[j][k] -= factor * A[i][k];
//            }
//            b[j] -= factor * b[i];
//        }
//    }
//
//    // Обратный ход метода Гаусса
//    for (int i = n - 1; i >= 0; --i) {
//        for (int j = i - 1; j >= 0; --j) {
//            b[j] -= A[j][i] * b[i];
//            A[j][i] = 0;
//        }
//    }
//}
//
//// Вычисление нормы матрицы
//double matrixNorm(const vector<vector<double>>& A) {
//    int n = A.size();
//    double maxNorm = 0;
//    for (int j = 0; j < n; ++j) {
//        double colNorm = 0;
//        for (int i = 0; i < n; ++i) {
//            colNorm += abs(A[i][j]);
//        }
//        maxNorm = max(maxNorm, colNorm);
//    }
//    return maxNorm;
//}
//
//// Вычисление обусловленности матрицы
//double conditionNumber(const vector<vector<double>>& A) {
//    double normA = matrixNorm(A);
//    vector<vector<double>> invA = inverseMatrix(A);
//    double normAInv = matrixNorm(invA);
//    return normA * normAInv;
//}
//
//// Вычисление невязки
//double residual(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
//    int n = A.size();
//    vector<double> Ax(n, 0);
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            Ax[i] += A[i][j] * x[j];
//        }
//    }
//    double res = 0;
//    for (int i = 0; i < n; ++i) {
//        res += abs(Ax[i] - b[i]);
//    }
//    return res;
//}
//
//int main() {
//    vector<vector<double>> A = { {1.00, 0.80, 0.64},
//                                {1.00, 0.90, 0.81},
//                                {1.00, 1.10, 1.21} };
//    vector<double> b = { erf(0.80), erf(0.90), erf(1.10) };
//
//    // Создаем копию матрицы A для вычисления обратной матрицы
//    vector<vector<double>> A_copy = A;
//
//    gaussElimination(A, b);
//
//    cout << "Solution: ";
//    for (double x : b) {
//        cout << x << " ";
//    }
//    cout << endl;
//
//    double cond = conditionNumber(A_copy);
//    cout << "Condition number: " << cond << endl;
//
//    double res = residual(A_copy, b, b);
//    cout << "Residual: " << res << endl;
//
//    double sumX = 0;
//    for (double x : b) {
//        sumX += x;
//    }
//    cout << "Sum of x: " << sumX << endl;
//    cout << "erf(1.0): " << erf(1.0) << endl;
//
//    vector<vector<double>> invA = inverseMatrix(A_copy);
//    cout << "Inverse matrix: " << endl;
//    for (const auto& row : invA) {
//        for (double val : row) {
//            cout << val << " ";
//        }
//        cout << endl;
//    }
//
//    return 0;
//}
