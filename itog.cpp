#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;
const double M_PI = 3.141592653589793;
const double L = M_PI / 2.0;

// Точное решение
double u_exact(double x, double t) {
    return exp(-t) * sin(x);
}

// Метод Прогонки для трёхдиагональной системы:
// Массивы индексируются с 1 до n;
void progonka( const vector<double>& a, const vector<double>& d, const vector<double>& c, const vector<double>& b, vector<double>& x) {
    int n = (int)d.size() - 1;
    vector<double> alpha(n + 1), beta(n + 1);

    // Прямой ход
    alpha[1] = d[1];
    beta[1] = b[1];
    for (int i = 2; i <= n; ++i) {
        double phi = a[i] / alpha[i - 1];
        alpha[i] = d[i] - phi * c[i - 1];
        beta[i] = b[i] - phi * beta[i - 1];
    }

    // Обратный ход
    x[n] = beta[n] / alpha[n];
    for (int i = n - 1; i >= 1; --i) {
        x[i] = (beta[i] - c[i] * x[i + 1]) / alpha[i];
    }
}

// Возвращает матрицу решений U для явной схемы
vector<vector<double>> Solexplicit(double T, int N, int M) {
    double h = L / N;
    double tau = T / M;
    double r = tau / (h * h);

    if (r > 0.5) {
        cerr << "Внимание: r = " << r << " > 0.5, явная схема неустойчива.\n";
    }

    vector<vector<double>> U(M + 1, vector<double>(N + 1, 0.0));

    // Начальный слой n=0:
    for (int i = 0; i <= N; ++i) {
        double x_i = i * h;
        U[0][i] = sin(x_i);
    }
    U[0][0] = 0.0;  // u(0,0) = 0
    U[0][N] = 1.0;  // u(L,0) = sin(L)=1, e^{-0}=1

    // Вспомогательный вектор для слоя n+1
    vector<double> u_next(N + 1, 0.0);

    for (int n = 0; n < M; ++n) {
        double t_np1 = (n + 1) * tau;

        // Граничные условия в слое n+1
        u_next[0] = 0.0;
        u_next[N] = exp(-t_np1);

        // Внутренние i=1..N-1
        for (int i = 1; i < N; ++i) {
            u_next[i] = U[n][i] + r * (U[n][i + 1] - 2.0 * U[n][i] + U[n][i - 1]);
        }

        // Копируем в U[n+1][...]
        for (int i = 0; i <= N; ++i) {
            U[n + 1][i] = u_next[i];
        }
    }

    return U;
}

// Возвращает матрицу решений U для неявной схемы
vector<vector<double>> Solimplicit(double T, int N, int M) {
    double h = L / N;
    double tau = T / M;
    double r = tau / (h * h);

    vector<vector<double>> U(M + 1, vector<double>(N + 1, 0.0));

    // Начальный слой
    for (int i = 0; i <= N; ++i) {
        double x_i = i * h;
        U[0][i] = sin(x_i);
    }
    U[0][0] = 0.0;
    U[0][N] = 1.0;

    int size = N - 1;  // число внутренних узлов
    vector<double> a(size + 1), d(size + 1), c(size + 1), b(size + 1), sol(size + 1);

    for (int n = 0; n < M; ++n) {
        double t_np1 = (n + 1) * tau;

        // Граничные
        U[n + 1][0] = 0.0;
        U[n + 1][N] = exp(-t_np1);

        // Заполняем a, d, c для всех внутренних i=1..N-1
        for (int i = 1; i <= size; ++i) {
            a[i] = -r;
            d[i] = 1.0 + 2.0 * r;
            c[i] = -r;
        }

        // Правая часть b[i]
        for (int i = 1; i <= size; ++i) {
            if (i == 1) {
                b[i] = U[n][i] + r * U[n + 1][0];  // U[n+1][0] = 0
            }
            else if (i == size) {
                b[i] = U[n][i] + r * U[n + 1][N];  // U[n+1][N] = e^{-t_np1}
            }
            else {
                b[i] = U[n][i];
            }
        }

        // Решаем тридиагональную систему методом Томаса
        progonka(a, d, c, b, sol);

        // Копируем sol[1..size] в U[n+1][1..N-1]
        for (int i = 1; i < N; ++i) {
            U[n + 1][i] = sol[i];
        }
    }
    return U;
}

void print_max_errors( const vector<vector<double>>& U, double T, int N, int M) {
    double h = L / N;
    double tau = T / M;
    int n0 = 0;
    int n1 = int(round(0.25 * M));
    int n2 = int(round(0.50 * M));
    int n3 = int(round(0.75 * M));
    int n4 = M;
    vector<int> layers = { n0, n1, n2, n3, n4 };
    vector<double> times;
    for (int n : layers) {
        times.push_back(n * tau);
    }

    cout << "\nМаксимальные ошибки на контрольных слоях:\n";
    cout << "  слой\tt^n\t\tmax|error|\n";
    for (int k = 0; k < (int)layers.size(); ++k) {
        int n = layers[k];
        double t_n = times[k];
        double max_err = 0.0;
        for (int i = 0; i <= N; ++i) {
            double x_i = i * h;
            double diff = fabs(U[n][i] - u_exact(x_i, t_n));
            if (diff > max_err) {
                max_err = diff;
            }
        }
        cout << setw(6) << n << "\t" << setw(6) << t_n;
        cout << "\t" << setw(12) << scientific << max_err << "\n";
    }
}

int main() {
    cout << fixed << setprecision(6);
    setlocale(LC_ALL, "RUS");
    // 1) Считываем T, N, M, схему
    double T;
    int N, M, type;
    cout << "Введите конечное время T: ";
    cin >> T;
    cout << "Введите число разбиений по x: N = ";
    cin >> N;
    cout << "Введите число шагов по t: M = ";
    cin >> M;
    cout << "Выберите схему: 1 - явная,  2 - неявная: ";
    cin >> type;

    // 2) Расчитаем h, tau, r
    double h = L / N;
    double tau = T / M;
    double r = tau / (h * h);
    cout << "\nПараметры сетки:\n";
    cout << "  L   = " << L << "\n";
    cout << "  h   = L/N  = " << h << "\n";
    cout << "  T   = " << T << "\n";
    cout << "  tau = T/M  = " << tau << "\n";
    cout << "  r   = tau/h^2 = " << r << "\n\n";

    // 3) В зависимости от выбора схемы строим матрицу U
    vector<vector<double>> U;
    if (type == 1) {
        cout << "Запускаем явную схемy...\n";
        U = Solexplicit(T, N, M);
    }
    else if (type == 2) {
        cout << "Запускаем неявную схему...\n";
        U = Solimplicit(T, N, M);
    }
    else {
        cerr << "Ошибка: выбор схемы должен быть 1 или 2.\n";
        return 1;
    }
    print_max_errors(U, T, N, M);

    // 4) Определяем индексы пяти слоёв: 0, 0.25T, 0.5T, 0.75T, T
    int n0 = 0;
    int n1 = int(round(0.25 * M));
    int n2 = int(round(0.50 * M));
    int n3 = int(round(0.75 * M));
    int n4 = M;

    // 5) Сохраняем их все в одном файле layers.txt
    ofstream fout("layers.txt");
    if (!fout.is_open()) {
        cerr << "Не удалось открыть файл layers.txt для записи\n";
        return 1;
    }

    fout << "x_i\tu_ex(0)\tu_num(0)\tu_ex(0.25T)\tu_num(0.25T)\tu_ex(0.50T)\tu_num(0.50T)\tu_ex(0.75T)\tu_num(0.75T)\tu_ex(T)\tu_num(T)\n";
    for (int i = 0; i <= N; ++i) {
        double x_i = i * h;
        double u_ex0 = u_exact(x_i, n0 * tau);
        double u_num0 = U[n0][i];

        double u_ex1 = u_exact(x_i, n1 * tau);
        double u_num1 = U[n1][i];

        double u_ex2 = u_exact(x_i, n2 * tau);
        double u_num2 = U[n2][i];

        double u_ex3 = u_exact(x_i, n3 * tau);
        double u_num3 = U[n3][i];

        double u_ex4 = u_exact(x_i, n4 * tau);
        double u_num4 = U[n4][i];

        fout << setw(8) << x_i << "\t"
            << setw(8) << u_ex0 << "\t" << setw(12) << u_num0 << "\t"
            << setw(8) << u_ex1 << "\t" << setw(12) << u_num1 << "\t"
            << setw(8) << u_ex2 << "\t" << setw(12) << u_num2 << "\t"
            << setw(8) << u_ex3 << "\t" << setw(12) << u_num3 << "\t"
            << setw(8) << u_ex4 << "\t" << setw(12) << u_num4 << "\n";
    }

    fout.close();
    cout << "\nВсе пять слоёв сохранены в файл layers.txt\n";
    return 0;
}
