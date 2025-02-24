//1
//#include <iostream>
//#include <cmath>
//#include <iomanip>
//#include <limits>
//
//using namespace std;
//
//const double PI = 3.14159265358979323846;
//
//double factorial(int n) {
//    if (n == 0 || n == 1) return 1;
//    return n * factorial(n - 1);
//}
//
//double My_erf(double x, int max_terms = 1000) {
//    double sum = 0.0;
//    double term = 0.0;
//    double prev_sum = 0.0;
//
//    for (int n = 0; n < max_terms; ++n) {
//        term = pow(-1, n) * pow(x, 2 * n + 1) / (factorial(n) * (2 * n + 1));
//        prev_sum = sum;
//        sum += term;
//
//        if (abs(sum - prev_sum) < numeric_limits<double>::epsilon()) {
//            break;
//        }
//    }
//
//    return (2.0 / sqrt(PI)) * sum;
//}
//
//int main() {
//    double x_values[] = { 0.5, 1.0, 5.0, 10.0 };
//
//    cout << fixed << setprecision(15);
//
//    for (int i = 0; i < 4; ++i) {
//        double x = x_values[i];
//        double calculated_erf = My_erf(x);
//        double library_erf = erf(x);
//
//        cout << "x = " << x << endl;
//        cout << "Calculated erf(x): " << calculated_erf << endl;
//        cout << "Library erf(x): " << library_erf << endl;
//        cout << "Difference: " << abs(calculated_erf - library_erf) << endl;
//        cout << "----------------------------------------" << endl;
//    }
//
//    return 0;
//}

//2

//#include <iostream>
//#include <cmath>
//#include<iomanip>
//using namespace std;
//
//double phi(double x, double epsilon) {
//    double result = 0;
//    double temp;
//    for (int k = 1;; k++) {
//        temp = 1 / (k * (k + x));
//        result += temp;
//        if (temp <= epsilon) break;
//    }
//
//    return result;
//}
//
//double difference(double x, double epsilon) {
//    double result = 0;
//    double temp;
//    for (int k = 1; ; k++) {
//        temp = (1 - x) / (k * (k + x) * (k + 1));
//        result += temp;
//        if (temp <= epsilon) break;
//    }
//    return result;
//}
//
//int main() {
//    setlocale(LC_ALL, "RUS");
//    cout << fixed << setprecision(5);
//    cout << "Доказательство, что phi(1)=1:" << endl;
//    double epsilon = 0.5e-8;
//    double step = 0.1;
//    double sum_phi = 0;
//    for (int k = 1; k <= 10000; k++) {
//        sum_phi += (1 / k) - (1 / (k + 1));
//    }
//    cout << "По формуле получим, что phi(1) = " << sum_phi << endl << "---------------------" << endl;
//
//    for (double x = 0.0; x <= 1.0; x += step) {
//        double result_phi = phi(x, epsilon);
//        double result_difference = difference(x, epsilon);
//
//        cout << "x = " << x << ", phi(x) = " << result_phi << ", difference(x) = " << result_difference << endl;
//    }
//
//    return 0;
//}

//3


//#include <iostream>
//#include <cmath>
//#include <iomanip>
//using namespace std;
//
//double testS(double x, double epsilon, int& n) {
//    double sum1 = 0, sum2 = 0;
//    for (int k = 1; ; k++) {
//        double temp1 = 1 / sqrt(pow(k, 3) + x);
//        double temp2 = 1 / sqrt(pow(k, 3) - x);
//        if (abs(temp1) < epsilon && abs(temp2) < epsilon) break;
//        sum1 += temp1;
//        sum2 += temp2;
//        n++;
//    }
//    double result = sum1 - sum2;
//    return result;
//}
//double S(double x, double epsilon, int& n) {
//    double sum = 0.0;
//    for (int k = 1; ; ++k) {
//        double temp = 1.0 / sqrt(k * k * k + x) - 1.0 / sqrt(k * k * k - x);
//        sum += temp;
//        n++;
//        if (abs(temp) < epsilon) break;
//    }
//    return sum;
//}
//
//double optimizedS(double x, double epsilon, int&n) {
//    double sum = 0.0;
//    for (int k = 1; ; ++k) {
//        double temp = -2.0 * x / (sqrt((k * k * k + x) * (k * k * k - x)) * (sqrt(k * k * k - x) + sqrt(k * k * k + x)));
//        sum += temp;
//        n++;
//        if (abs(temp) < epsilon) break;
//    }
//    return sum;
//}
//
//int main() {
//    cout << fixed << setprecision(6);
//    setlocale(LC_ALL, "RUS");
//    double epsilon = 3e-8;
//    int n1 = 0, n2 = 0, on1 = 0, on2 = 0;
//    double x1 = 0.5, x2 = 0.999999999;
//    double f1, f2, of1, of2;
//    f1 = testS(x1, epsilon, n1);
//    f2 = testS(x2, epsilon, n2);
//    of1 = optimizedS(x1, epsilon, on1);
//    of2 = optimizedS(x2, epsilon, on2);
//    cout << "Исходная функция:" << endl;
//    cout << "x=0.5" << endl << "Количество итераций: " << n1 << endl << "Время: " << n1 * 500 << " mc" << endl<<"Результат: "<<abs(f1)<<endl;
//    cout << "--------------------------" << endl;
//    cout<<"x=0.999999999"<<endl<< "Количество итераций: " << n2 << endl << "Время: " << n2 * 500 << " mc" << endl << "Результат: " << abs(f2) << endl;
//    cout << "--------------------------" << endl;
//    cout << "Оптимизированная функция:" << endl;
//    cout << "x=0.5" << endl << "Количество итераций: " << on1 << endl << "Время: " << on1 * 500 << " mc" << endl << "Результат: " << abs(of1) << endl;
//    cout << "--------------------------" << endl;
//    cout << "x=0.999999999" << endl << "Количество итераций: " << on2 << endl << "Время: " << on2 * 500 << " mc" << endl << "Результат: " << abs(of2) << endl;
//    cout << "--------------------------" << endl;
//    return 0;
//}

//4

//#include <iostream>
//#include <cmath>
//#include <iomanip>
//#define _USE_MATH_DEFINES
//#include <math.h>
//
//using namespace std;
//
//double sum_first(double epsilon, long long& iter) {
//    double sum = 0;
//    for (int n = 1; ; ++n) {
//        double temp = 1.0 / (n * n + 1.0);
//        sum += temp;
//        iter++;
//        if (temp < epsilon) break;
//    }
//
//    return sum;
//}
//
//double sum_second(double epsilon, long long& iter) {
//    double sum = 0;
//    for (int n = 1; ; ++n) {
//        double temp = 1.0 / (pow(n, 4) * (n * n + 1.0));
//        sum += temp ;
//        iter++;
//        if (temp <= epsilon) break; 
//    }
//
//    return pow(M_PI, 2) / 6 - pow(M_PI, 4) / 90 + sum;
//}
//int main() {
//    setlocale(LC_ALL, "RUS");
//    double epsilon = 1e-10;
//    long long iterations1 = 0, iterations2 = 0;
//    double sum1 = 0, sum2 = 0;
//    sum1 = sum_first(epsilon, iterations1);
//    sum2 = sum_second(epsilon, iterations2);
//    cout.precision(6);
//    cout << "Первый способ: " << endl;
//    cout << "Сумма ряда: " << sum1 << endl;
//    cout << "Количество итераций: " << iterations1 << endl; 
//    cout << "-------------------------" << endl;
//    cout << "Второй способ: " << endl;
//    cout << "Сумма ряда: " << sum2 << endl;
//    cout << "Количество итераций: " << iterations2 << endl;
//    return 0;
//}
