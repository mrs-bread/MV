#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <iomanip>

using namespace std;
const double M_PI = acos(-1.0);
const double M_PI_2 = M_PI / 2.0;
mt19937_64 rng(random_device{}());
const long N = 100;
void Графики() {
    ofstream f1("Задача1.txt");
    for (double x = 0; x <= 20.0; x += 0.1) {
        double y = (x < 3.0 ? 10.0 * x / 3.0 : 10.0 * (20.0 - x) / 17.0);
        f1 << x << " " << y << "\n";
    }
    f1.close();

    ofstream f2("Задача2.txt");
    for (double x = 0; x <= 5.0; x += 0.01) {
        double y = sqrt(11.0 - 3.0 * pow(sin(x), 2));
        f2 << x << " " << y << "\n";
    }
    f2.close();

    ofstream f3("Задача3.txt");
    double R = 3.0;
    for (double phi = 0; phi <= 2 * M_PI; phi += 0.01) {
        double xx = R * cos(phi);
        double yy = R * sin(phi);
        f3 << xx << " " << yy << "\n";
    }
    f3.close();

    ofstream f4("Задача4.txt");
    for (double phi = 0; phi <= 2 * M_PI; phi += 0.01) {
        double rho = sqrt(14.0 * cos(phi) * cos(phi)
            + 8.0 * sin(phi) * sin(phi));
        double xx = rho * cos(phi);
        double yy = rho * sin(phi);
        f4 << xx << " " << yy << "\n";
    }
    f4.close();
}

pair<double, long> M_K_1() {
    ofstream pts("Задача1_Точки.txt");
    uniform_real_distribution<double> dx(0.0, 20.0), dy(0.0, 10.0);
    long inside = 0;
    for (long i = 0; i < N; ++i) {
        double x = dx(rng), y = dy(rng);
        double fx = (x < 3.0 ? 10.0 * x / 3.0 : 10.0 * (20.0 - x) / 17.0);
        bool in = (y <= fx);
        pts << x << " " << y << " " << (in ? 1 : 0) << "\n";
        if (in) ++inside;
    }
    pts.close();
    double areaR = 20.0 * 10.0;
    return { areaR * inside / double(N), inside };
}

pair<double, long> M_K_2() {
    ofstream pts("Задача2_Точки.txt");
    uniform_real_distribution<double> dx(0.0, 5.0), dy(0.0, sqrt(11.0));
    long inside = 0;
    for (long i = 0; i < N; ++i) {
        double x = dx(rng), y = dy(rng);
        bool in = (y <= sqrt(11.0 - 3.0 * pow(sin(x), 2)));
        pts << x << " " << y << " " << (in ? 1 : 0) << "\n";
        if (in) ++inside;
    }
    pts.close();
    double areaR = 5.0 * sqrt(11.0);
    return { areaR * inside / double(N), inside };
}

pair<double, long> M_K_3() {
    ofstream pts("Задача3_Точки.txt");
    double R = 3.0;
    uniform_real_distribution<double> dxy(-R, R);
    long inside = 0;
    for (long i = 0; i < N; ++i) {
        double x = dxy(rng), y = dxy(rng);
        bool in = (x * x + y * y <= R * R);
        pts << x << " " << y << " " << (in ? 1 : 0) << "\n";
        if (in) ++inside;
    }
    pts.close();
    return { 4.0 * inside / double(N), inside };
}

pair<double, long> M_K_4() {
    ofstream pts("Задача4_Точки.txt");
    double A = 14.0, B = 8.0;
    double xmax = sqrt(A), ymax = sqrt(B);
    uniform_real_distribution<double> dx(-xmax, xmax), dy(-ymax, ymax);
    long inside = 0;
    for (long i = 0; i < N; ++i) {
        double x = dx(rng), y = dy(rng);
        double phi = atan2(y, x);
        bool in = (hypot(x, y) <= sqrt(A * cos(phi) * cos(phi) + B * sin(phi) * sin(phi)));
        pts << x << " " << y << " " << (in ? 1 : 0) << "\n";
        if (in) ++inside;
    }
    pts.close();
    double areaR = 2 * xmax * 2 * ymax;
    return { areaR * inside / double(N), inside };
}

int main() {
    cout << fixed << setprecision(6);
    setlocale(LC_ALL, "RUS");
    Графики();
    double true1 = 15.0 + 85.0;
    pair<double, long> result1 = M_K_1();
    double area1 = result1.first;
    long   попало1 = result1.second;
    double err1 = fabs(area1 - true1);
    cout << "Задача 1:\n"
        << "  Попало = " << setw(7) << попало1 << " из " << N << "\n"
        << "  Оценка площади = " << setw(9) << area1 << "\n"
        << "  Аналитическое = " << setw(7) << true1 << "\n"
        << "  Абсолютная погрешность = " << setw(8) << err1 << "\n"
        << "  Относительная погрешность = " << setw(6) << (err1 / true1 * 100) << "%\n\n";

    double true2 = 15.3194;
    pair<double, long> result2 = M_K_2();
    double I2 = result2.first;
    long   попало2 = result2.second;
    double err2 = fabs(I2 - true2);
    cout << "Задача 2:\n"
        << "  Попало = " << setw(7) << попало2 << " из " << N << "\n"
        << "  Оценка интеграла = " << setw(9) << I2 << "\n"
        << "  Аналитическое = " << setw(7) << true2 << "\n"
        << "  Абсолютная погрешность = " << setw(8) << err2 << "\n"
        << "  Относительная погрешность = " << setw(6) << (err2 / true2 * 100) << "%\n\n";

    pair<double, long> result3 = M_K_3();
    double pi_estimate = result3.first;
    long   попало3 = result3.second;
    double err3 = fabs(pi_estimate - M_PI);
    cout << "Задача 3:\n"
        << "  Попало = " << setw(7) << попало3 << " из " << N << "\n"
        << "  Оценка Пи   = " << setw(9) << pi_estimate << "\n"
        << "  Аналитическое = " << setw(8) << M_PI << "\n"
        << "  Абсолютная погрешность = " << setw(8) << err3 << "\n"
        << "  Относительная погрешность = " << setw(6) << (err3 / M_PI * 100) << "%\n\n";

    double true4 = M_PI * sqrt(14.0 * 8.0);
    pair<double, long> result4 = M_K_4();
    double area4 = result4.first;
    long   попало4 = result4.second;
    double err4 = fabs(area4 - true4);
    cout << "Задача 4:\n"
        << "  Попало = " << setw(7) << попало4 << " из " << N << "\n"
        << "  Оценка площади = " << setw(9) << area4 << "\n"
        << "  Аналитическое = " << setw(7) << true4 << "\n"
        << "  Абсолютная погрешность = " << setw(8) << err4 << "\n"
        << "  Относительная погрешность = " << setw(6) << (err4 / true4 * 100) << "%\n";

    return 0;
}
