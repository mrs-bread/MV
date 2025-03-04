#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Функция для вывода матрицы
void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
}

// Метод Гаусса для приведения матрицы к ступенчатому виду
int gaussianElimination(vector<vector<double>>& matrix) {
    int rank = 0;
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int col = 0, row = 0; col < cols && row < rows; ++col) {
        // Находим максимальный элемент в текущем столбце
        int pivotRow = row;
        for (int i = row + 1; i < rows; ++i) {
            if (fabs(matrix[i][col]) > fabs(matrix[pivotRow][col])) {
                pivotRow = i;
            }
        }

        // Если все элементы в столбце нулевые, переходим к следующему столбцу
        if (fabs(matrix[pivotRow][col]) < 1e-9) {
            continue;
        }

        // Меняем строки местами
        swap(matrix[row], matrix[pivotRow]);

        // Обнуляем элементы ниже главной диагонали
        for (int i = row + 1; i < rows; ++i) {
            double factor = matrix[i][col] / matrix[row][col];
            for (int j = col; j < cols; ++j) {
                matrix[i][j] -= factor * matrix[row][j];
            }
        }

        ++row;
        ++rank;
    }

    return rank;
}

int main() {
    // Определяем матрицу A и вектор b
    vector<vector<double>> A = {
        {0.1, 0.2, 0.3},
        {0.4, 0.5, 0.6},
        {0.7, 0.8, 0.9}
    };

    vector<double> b = {0.1, 0.3, 0.5};

    // Создаем расширенную матрицу (A|b)
    vector<vector<double>> Ab(A.size(), vector<double>(A[0].size() + 1));
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) {
            Ab[i][j] = A[i][j];
        }
        Ab[i][A[0].size()] = b[i];
    }

    cout << "Исходная расширенная матрица (A|b):" << endl;
    printMatrix(Ab);

    // Приводим матрицы к ступенчатому виду и вычисляем ранги
    vector<vector<double>> A_copy = A; // Копия матрицы A
    int rankA = gaussianElimination(A_copy);
    int rankAb = gaussianElimination(Ab);

    cout << "\nRank(A): " << rankA << endl;
    cout << "Rank(A|b): " << rankAb << endl;

    // Анализ результатов
    int n = A[0].size(); // количество неизвестных
    if (rankA == rankAb && rankA < n) {
        cout << "Система имеет бесконечное множество решений." << endl;
    } else if (rankA == rankAb && rankA == n) {
        cout << "Система имеет единственное решение." << endl;
    } else {
        cout << "Система несовместна (решений нет)." << endl;
    }

    return 0;
}
