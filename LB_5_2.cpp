#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
using namespace std;
pair<double, double> Крамер2x2(double a11, double a12, double a21, double a22, double b1, double b2) {
    double det = a11 * a22 - a12 * a21;
    if (fabs(det) < 1e-12) {
        cerr << "Вырожденная матрица" << endl;
        return make_pair(0.0, 0.0);
    }
    double u = (b1 * a22 - a12 * b2) / det;
    double v = (a11 * b2 - b1 * a21) / det;
    return make_pair(u, v);
}

vector<double> Крамер3x3(const vector<vector<double>>& A, const vector<double>& B) {
    double a11 = A[0][0], a12 = A[0][1], a13 = A[0][2];
    double a21 = A[1][0], a22 = A[1][1], a23 = A[1][2];
    double a31 = A[2][0], a32 = A[2][1], a33 = A[2][2];
    double det = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
    if (fabs(det) < 1e-12) {
        cerr << "Вырожденная матрица" << endl;
        return vector<double>(3, 0.0);
    }
    double d1 = B[0] * (a22 * a33 - a23 * a32) - a12 * (B[1] * a33 - a23 * B[2]) + a13 * (B[1] * a32 - a22 * B[2]);
    double d2 = a11 * (B[1] * a33 - a23 * B[2]) - B[0] * (a21 * a33 - a23 * a31) + a13 * (a21 * B[2] - B[1] * a31);
    double d3 = a11 * (a22 * B[2] - B[1] * a32) - a12 * (a21 * B[2] - B[1] * a31) + B[0] * (a21 * a32 - a22 * a31);
    vector<double> sol(3);
    sol[0] = d1 / det;
    sol[1] = d2 / det;
    sol[2] = d3 / det;
    return sol;
}

int main() {
    setlocale(LC_ALL, "RUS");
    int n = 6;
    vector<double> x = { 3,5,7,9,11,13 };
    vector<double> y = { 26,76,150,240,360,500 };
    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    double Slnx = 0, Slny = 0, Slnln = 0, Slnylnx = 0, Sxlny = 0;
    double Sx2 = 0, Sx3 = 0, Sx4 = 0, Sx2y = 0;
    for (int i = 0; i < n; ++i) {
        Sx += x[i];
        Sy += y[i];
        Sxx += x[i] * x[i];
        Sxy += x[i] * y[i];
        double lx = log(x[i]);
        double ly = log(y[i]);
        Slnx += lx;
        Slny += ly;
        Slnln += lx * lx;
        Slnylnx += lx * ly;
        Sxlny += x[i] * ly;
        double x2 = x[i] * x[i];
        Sx2 += x2;
        Sx3 += x2 * x[i];
        Sx4 += x2 * x2;
        Sx2y += x2 * y[i];
    }
    pair<double, double> linf = Крамер2x2(Sxx, Sx, Sx, n, Sxy, Sy);
    double a_lin = linf.first;
    double b_lin = linf.second;
    pair<double, double> powf = Крамер2x2(Slnln, Slnx, Slnx, n, Slnylnx, Slny);
    double beta_pow = powf.first;
    double a_pow = exp(powf.second);
    pair<double, double> expf = Крамер2x2(Sxx, Sx, Sx, n, Sxlny, Slny);
    double beta_exp = expf.first;
    double a_exp = exp(expf.second);
    vector<vector<double>> A = { {Sx4, Sx3, Sx2},
                                {Sx3, Sx2, Sx },
                                {Sx2, Sx , (double)n} };
    vector<double> B = { Sx2y, Sxy, Sy };
    vector<double> quadf = Крамер3x3(A, B);
    double a_quad = quadf[0];
    double b_quad = quadf[1];
    double c_quad = quadf[2];

    auto r2 = [](double v) { return round(v * 100.0) / 100.0; };
    a_lin = r2(a_lin);   b_lin = r2(b_lin);
    a_pow = r2(a_pow);   beta_pow = r2(beta_pow);
    a_exp = r2(a_exp);   beta_exp = r2(beta_exp);
    a_quad = r2(a_quad);  b_quad = r2(b_quad);
    c_quad = r2(c_quad);

    cout << fixed << setprecision(2);
    cout << "Линейная:         y = " << a_lin << " * x + " << b_lin << endl;
    cout << "Степенная:        y = " << a_pow << " * x^" << beta_pow << endl;
    cout << "Экспоненциальная: y = " << a_exp << " * e^(" << beta_exp << " * x)" << endl;
    cout << "Квадратичная:     y = " << a_quad << " * x^2 + " << b_quad << " * x + " << c_quad << endl;
    cout <<"----------------------------------------------------------------" << endl;

    double E_lin = 0, E_pow = 0, E_exp = 0, E_quad = 0;
    for (int i = 0; i < n; ++i) {
        double y_lin = a_lin * x[i] + b_lin;
        double y_pow = a_pow * pow(x[i], beta_pow);
        double y_exp = a_exp * exp(beta_exp * x[i]);
        double y_quad = a_quad * x[i] * x[i] + b_quad * x[i] + c_quad;
        E_lin += pow(y[i] - y_lin, 2);
        E_pow += pow(y[i] - y_pow, 2);
        E_exp += pow(y[i] - y_exp, 2);
        E_quad += pow(y[i] - y_quad, 2);
    }
    cout << "Суммарные погрешности:" << endl;
    cout << fixed << setprecision(4);
    cout << "Линейная = " << E_lin << endl 
        << "Степенная = " << E_pow << endl 
        << "Экспоненциальная = " << E_exp << endl 
        << "Квадратичная = " << E_quad << endl;
    cout << "----------------------------------------------------------------" << endl;
    double minS = min({ E_lin, E_pow, E_exp, E_quad });
    string best;
    if (minS == E_lin)  best = "Линейная";
    else if (minS == E_pow)  best = "Степенная";
    else if (minS == E_exp)  best = "Экспоненциальная";
    else best = "Квадратичная";
    cout << "Лучшая апроксимирующая функция: " << best << " функция."<<endl;
    cout << "----------------------------------------------------------------" << endl;
    ofstream Result("all_data.txt");
    Result << "x\tИсходная\tЛинейная\tСтепенная\tЭкспоненциальная\tКвадратичная\n";
    for (int i = 0; i < n; ++i) {
        double y_lin = a_lin * x[i] + b_lin;
        double y_pow = a_pow * pow(x[i], beta_pow);
        double y_exp = a_exp * exp(beta_exp * x[i]);
        double y_quad = a_quad * x[i] * x[i] + b_quad * x[i] + c_quad;
        Result << x[i] << '\t' << y[i] << '\t'
            << y_lin << '\t' << y_pow << '\t'
            << y_exp << '\t' << y_quad << '\n';
    }
    Result.close();
    cout << "Данные сохранены." << endl;
    return 0;
}
