#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;
//Если для каждой строки i выполняется |A[i][i]| > ∑(|A[i][j]|, j != i)
bool ДиагональноеПреобладание(const vector<vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double diag = fabs(A[i][i]);
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                sum += fabs(A[i][j]);
            }
        }
        if (diag < sum) {
            return false;
        }
    }
    return true;
}

double НормаНевязки(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
    int n = b.size();
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        double diff = sum - b[i];
        result += pow(diff,2);
    }
    return sqrt(result);
}

// x_[i] = (b[i] - ∑(j!=i) A[i][j] * x[j]) / A[i][i]
vector<double> Якоби(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& x0, double Точность, int МаксИтераций, const string& filename)
{
    int n = b.size();
    vector<double> x = x0;
    vector<double> x_new(n, 0.0);

    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Не удалось открыть файл " << filename << endl;
        return x;
    }

    int iter = 0;
    double result = НормаНевязки(A, x, b);
    outFile << iter << " " << result << "\n";

    while (result > Точность && iter < МаксИтераций) {
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sigma += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sigma) / A[i][i];
        }

        x = x_new;
        iter++;

        result = НормаНевязки(A, x, b);
        outFile << iter << " " << result << "\n";
    }

    outFile.close();
    cout << "Метод Якоби: достигнута точность " << Точность << " на итерации " << iter << endl << "Приблизительное решение: ";
    for (double xi : x) {
        cout << xi << " ";
    }
    cout << "\nКонечная норма невязки: " << result << "\n" << endl;
    return x;
}

// x[i] = ( b[i] - ∑(j<i) A[i][j]*x[j] - ∑(j>i) A[i][j]*x[j] ) / A[i][i]
vector<double> Зейдель(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& x0, double Точность, int МаксИтераций, const string& filename)
{
    int n = b.size();
    vector<double> x = x0;
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Не удалось открыть файл " << filename << endl;
        return x;
    }
    int iter = 0;
    double result = НормаНевязки(A, x, b);
    outFile << iter << " " << result << "\n";

    while (result > Точность && iter < МаксИтераций) {
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < i; j++) {
                sigma += A[i][j] * x[j];
            }
            for (int j = i + 1; j < n; j++) {
                sigma += A[i][j] * x[j];
            }
            x[i] = (b[i] - sigma) / A[i][i];
        }
        iter++;
        result = НормаНевязки(A, x, b);
        outFile << iter << " " << result << "\n";
    }

    outFile.close();
    cout << "Метод Зейделя: достигнута точность " << Точность << " на итерации " << iter << endl << "Приблизительное решение: ";
    for (double xi : x) {
        cout << xi << " ";
    }
    cout << "\nКонечная норма невязки: " << result << "\n" << endl;

    return x;
}

int main() {
    setlocale(LC_ALL, "RUS");
    vector<vector<double>> A = {
        {12.14,  1.32,  -0.78, -2.75},
        {-0.89, 16.57,   1.88, -1.55},
        { 2.65, -1.27, -15.64, -0.64},
        { 2.44,  1.52,   1.93,-11.43}
    };
    vector<double> b = { 14.78, -12.14, -11.65, 4.26 };

    if (!ДиагональноеПреобладание(A)) {
        cout << "Матрица A не удовлетворяет условию диагонального преобладания. Итерационные методы могут не сходиться." << endl;
    }
    else {
        cout << "Матрица удовлетворяет условию диагонального преобладания." << endl;
    }
    double Точность = 1e-10;
    int МаксИтераций = 1000;

    cout << "\n=== Решение начального приближения (0,0,0,0) ===" << endl;
    vector<double> x0 = { 0.0, 0.0, 0.0, 0.0 };
    vector<double> Якоби_0 = Якоби(A, b, x0, Точность, МаксИтераций, "Якоби_0.txt");
    vector<double> Зейдель_0 = Зейдель(A, b, x0, Точность, МаксИтераций, "Зейдель_0.txt");

    cout << "\n=== Решение начального приближения (1,1,1,1) ===" << endl;
    vector<double> x0_1 = { 1.0, 1.0, 1.0, 1.0 };
    vector<double> Якоби_1 = Якоби(A, b, x0_1, Точность, МаксИтераций, "Якоби_1.txt");
    vector<double> Зейдель_1 = Зейдель(A, b, x0_1, Точность, МаксИтераций, "Зейдель_1.txt");

    cout << "\n==== Решение начального приближения (2,2,2,2) ====" << endl;
    vector <double> x0_2 = { 2.0,2.0,2.0,2.0 };
    vector<double> Якоби_2 = Якоби(A, b, x0_2, Точность, МаксИтераций, "Якоби_2.txt");
    vector<double> Зейдель_2 = Зейдель(A, b, x0_2, Точность, МаксИтераций, "Зейдель_2.txt");
    return 0;
}
