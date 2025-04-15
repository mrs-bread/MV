#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;
double Pi = 3.1415926535;
double f(double x) {
    return 1.0 / (1 + 25 * x * x);
}

double lagrangeInterpolation(double x, const vector<double>& nodes, const vector<double>& values) {
    double result = 0.0;
    int n = nodes.size();
    for (int i = 0; i < n; i++) {
        double term = values[i];
        for (int j = 0; j < n; j++) {
            if (j != i) {
                term *= (x - nodes[j]) / (nodes[i] - nodes[j]);
            }
        }
        result += term;
    }
    return result;
}

vector<double> equallySpacedNodes(double a, double b, int n) {
    vector<double> nodes(n);
    for (int i = 0; i < n; i++) {
        nodes[i] = a + i * (b - a) / (n - 1);
    }
    return nodes;
}

// Формула: x_i = (a+b)/2 + (b-a)/2 * cos((2i+1)*pi/(2n))
vector<double> chebyshevNodes(double a, double b, int n) {
    vector<double> nodes(n);
    for (int i = 0; i < n; i++) {
        double xi = cos((2 * i + 1) * Pi / (2 * n));
        nodes[i] = 0.5 * (a + b) + 0.5 * (b - a) * xi;
    }
    return nodes;
}

// Структура для хранения коэффициентов сегмента кубического сплайна
struct SplineSegment {
    double a, b, c, d; // коэффициенты полинома: a + b*(x-x_i) + c*(x-x_i)^2 + d*(x-x_i)^3
    double x;          // x-координата начала сегмента
};

// Вычисление коэффициентов натурального кубического сплайна по таблице (x, y)
vector<SplineSegment> computeCubicSpline(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    vector<double> alpha(n, 0.0);
    for (int i = 1; i < n - 1; i++) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    vector<double> l(n), mu(n), z(n);
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    for (int i = 1; i < n - 1; i++) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n - 1] = 1.0;
    z[n - 1] = 0.0;

    vector<double> c(n, 0.0), b(n - 1), d(n - 1);
    for (int j = n - 2; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    vector<SplineSegment> spline(n - 1);
    for (int i = 0; i < n - 1; i++) {
        spline[i].a = y[i];
        spline[i].b = b[i];
        spline[i].c = c[i];
        spline[i].d = d[i];
        spline[i].x = x[i];
    }
    return spline;
}

// Функция для вычисления значения кубического сплайна в заданной точке x_val
double evaluateSpline(const vector<SplineSegment>& spline, double x_val) {
    int n = spline.size();
    SplineSegment seg;
    if (x_val <= spline[0].x) {
        seg = spline[0];
    }
    else if (x_val >= spline[n - 1].x) {
        seg = spline[n - 1];
    }
    else {
        int i = 0;
        while (i < n - 1 && x_val >= spline[i + 1].x)
            i++;
        seg = spline[i];
    }
    double dx = x_val - seg.x;
    return seg.a + seg.b * dx + seg.c * dx * dx + seg.d * dx * dx * dx;
}

int main() {
    // ---------------------------
    // Часть 1. Интерполяция многочленом Лагранжа для функции Рунге на отрезке [-1, 1]
    // ---------------------------
    double a = -1.0, b = 1.0;
    int n_nodes = 11;  // количество узлов
    int n_points = 201; // число точек для построения графика

    // Равноотстоящие узлы
    vector<double> nodes_equal = equallySpacedNodes(a, b, n_nodes);
    vector<double> values_equal(n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        values_equal[i] = f(nodes_equal[i]);
    }

    // Чебышевские узлы
    vector<double> nodes_cheb = chebyshevNodes(a, b, n_nodes);
    vector<double> values_cheb(n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        values_cheb[i] = f(nodes_cheb[i]);
    }

    ofstream out_interp("interpolation_results.txt");
    out_interp << "x\tf(x)\tLagrange_equal\tLagrange_cheb\n";
    for (int i = 0; i < n_points; i++) {
        double x = a + i * (b - a) / (n_points - 1);
        double fx = f(x);
        double L_equal = lagrangeInterpolation(x, nodes_equal, values_equal);
        double L_cheb = lagrangeInterpolation(x, nodes_cheb, values_cheb);
        out_interp << setprecision(8) << x << "\t" << fx << "\t" << L_equal << "\t" << L_cheb << "\n";
    }
    out_interp.close();

    // ---------------------------
    // Часть 2. Кубический сплайн для функции Рунге на отрезке [-1, 1]
    // ---------------------------
    vector<double> nodes_spline = equallySpacedNodes(a, b, n_nodes);
    vector<double> values_spline(n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        values_spline[i] = f(nodes_spline[i]);
    }
    vector<SplineSegment> spline = computeCubicSpline(nodes_spline, values_spline);

    ofstream out_spline("spline_results.txt");
    out_spline << "x\tf(x)\tSpline\n";
    for (int i = 0; i < n_points; i++) {
        double x = a + i * (b - a) / (n_points - 1);
        double fx = f(x);
        double Sx = evaluateSpline(spline, x);
        out_spline << setprecision(8) << x << "\t" << fx << "\t" << Sx << "\n";
    }
    out_spline.close();

    // ---------------------------
    // Часть 3. Кубический сплайн для данных, заданных таблицей
    // Табличные данные: x: 2, 3, 5, 7 и f(x): 4, -2, 6, -3
    // ---------------------------
    vector<double> x_table = { 2.0, 3.0, 5.0, 7.0 };
    vector<double> y_table = { 4.0, -2.0, 6.0, -3.0 };
    vector<SplineSegment> spline_table = computeCubicSpline(x_table, y_table);

    // Для построения графика создаём плотную сетку от 2 до 7
    int n_table_points = 201;
    ofstream out_spline_table("spline_table_results.txt");
    out_spline_table << "x\tSpline\tOriginal\n";
    for (int i = 0; i < n_table_points; i++) {
        double x = x_table.front() + i * (x_table.back() - x_table.front()) / (n_table_points - 1);
        double Sx = evaluateSpline(spline_table, x);
        // Определяем, является ли текущая точка узловой (с заданной точностью)
        double orig = NAN;
        for (double xt : x_table) {
            if (fabs(x - xt) < 1e-8) {
                orig = xt; // здесь можно выводить исходное значение y_table, но для простоты напишем NaN, так как графически узловые точки можно отметить иначе
            }
        }
        // Записываем x, значение сплайна и (для контроля) значение исходных данных, если x совпадает с узлом
        // Если точка не узловая, в столбце Original будет NAN (не число).
        int idx = -1;
        for (size_t j = 0; j < x_table.size(); j++) {
            if (fabs(x - x_table[j]) < 1e-8) { idx = j; break; }
        }
        if (idx != -1)
            out_spline_table << setprecision(8) << x << "\t" << Sx << "\t" << y_table[idx] << "\n";
        else
            out_spline_table << setprecision(8) << x << "\t" << Sx << "\t" << "NaN" << "\n";
    }
    out_spline_table.close();

    cout << "Результаты сохранены в файлы:" << endl;
    cout << " - interpolation_results.txt" << endl;
    cout << " - spline_results.txt" << endl;
    cout << " - spline_table_results.txt" << endl;

    return 0;
}
