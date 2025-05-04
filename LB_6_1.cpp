// erf_table.cpp
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
const double M_PI = 3.1415926535;
// Интегрируем функцию exp(-t*t) на [0, x] по правилу Симпсона
double erf_simpson(double x, int n = 1000) {
    if (n % 2 != 0) n++;            // n должно быть чётным
    double h = x / n;
    double sum = exp(-0 * 0) + exp(-x * x);
    for (int i = 1; i < n; i++) {
        double t = i * h;
        sum += (i % 2 == 0 ? 2 : 4) * exp(-t * t);
    }
    return (2.0 / sqrt(M_PI)) * (h / 3.0) * sum;
}

int main() {
    cout << "   x\t\terf(x)\n";
    cout << "-----------------------\n";
    for (int k = 0; k <= 20; k++) {
        double x = k * 0.1;
        double val = erf_simpson(x);
        cout << fixed << setprecision(1) << x << "\t"
            << setprecision(6) << val << "\n";
    }
    return 0;
}


// //pi_approx.cpp
//#include <iostream>
//#include <cmath>
//#include <iomanip>
//using namespace std;
//
//double f(double x) { return 4.0 / (1 + x * x); }
//
//double rect_left(int n) {
//    double h = 1.0 / n, sum = 0;
//    for (int i = 0; i < n; i++) sum += f(i * h);
//    return h * sum;
//}
//
//double trapezoid(int n) {
//    double h = 1.0 / n, sum = 0.5 * (f(0) + f(1));
//    for (int i = 1; i < n; i++) sum += f(i * h);
//    return h * sum;
//}
//
//double simpson(int n) {
//    if (n % 2 != 0) return simpson(n + 1);
//    double h = 1.0 / n, sum = f(0) + f(1);
//    for (int i = 1; i < n; i++) {
//        sum += (i % 2 ? 4 : 2) * f(i * h);
//    }
//    return sum * h / 3.0;
//}
//
//int main() {
//    int ns[] = { 8, 32, 128 };
//    cout << " n\tRectangle\tTrapezoid\tSimpson\n";
//    cout << "---------------------------------------------\n";
//    for (int n : ns) {
//        cout << n << "\t"
//            << fixed << setprecision(10)
//            << rect_left(n) << "\t"
//            << trapezoid(n) << "\t"
//            << simpson(n) << "\n";
//    }
//    return 0;
//}


//// integral_strategy.cpp
//#include <iostream>
//#include <cmath>
//#include <iomanip>
//using namespace std;
//
//const double M_PI = 3.1415926535;
//// Метод Симпсона для func на [a,b] с n панели (n чётное)
//double simpson(double (*func)(double), double a, double b, int n) {
//    if (n % 2 != 0) ++n;
//    double h = (b - a) / n;
//    double sum = func(a) + func(b);
//    for (int i = 1; i < n; ++i) {
//        double x = a + i * h;
//        sum += (i % 2 == 0 ? 2.0 : 4.0) * func(x);
//    }
//    return sum * h / 3.0;
//}
//
//// Метод трапеций для func на [a,b] с n панели
//double trapezoid(double (*func)(double), double a, double b, int n) {
//    double h = (b - a) / n;
//    double sum = 0.5 * (func(a) + func(b));
//    for (int i = 1; i < n; ++i) {
//        double x = a + i * h;
//        sum += func(x);
//    }
//    return sum * h;
//}
//
//// f1(x) = e^{x^2}
//double f1(double x) {
//    return exp(x * x);
//}
//
//// f2(x) = 1 / (4 - sin(16π x))
//double f2(double x) {
//    return 1.0 / (4.0 - sin(16.0 * M_PI * x));
//}
//
//int main() {
//    // Число панелей
//    int n1 = 200;    // для Симпсона на [0,2]
//    int n2 = 400;    // для трапеций на [2,4] (лучше кратно периоду: 2 периода * 200 точек)
//
//    // Интегралы по частям
//    double I1 = simpson(f1, 0.0, 2.0, n1);
//    double I2 = trapezoid(f2, 2.0, 4.0, n2);
//    double I = I1 + I2;
//
//    cout << fixed << setprecision(10);
//    cout << "I1 = ∫0^2 e^{x^2} dx by Simpson  = " << I1 << "\n";
//    cout << "I2 = ∫2^4 1/(4 - sin(16πx)) dx by Trapezoid = " << I2 << "\n";
//    cout << "I  = I1 + I2                         = " << I << "\n";
//
//    return 0;
//}
