#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Точное решение
double exactSolution(double x) {
    return exp(sin(x));
}

// Правая часть ОДУ y' = y * cos(x)
double f(double y, double x) {
    return y * cos(x);
}

// Разностная схема первого порядка
vector<double> firstOrderScheme(double a, double b, double h, double y0) {
    int n = (b - a) / h + 1;
    vector<double> y(n);
    y[0] = y0;
    for (int i = 1; i < n; ++i) {
        double x = a + (i - 1) * h;
        y[i] = y[i - 1] + h * f(y[i - 1], x);
    }
    return y;
}

// Симметричная схема второго порядка
vector<double> secondOrderScheme(double a, double b, double h, double y0) {
    int n = (b - a) / h + 1;
    vector<double> y(n);
    y[0] = y0;
    for (int i = 1; i < n; ++i) {
        double x = a + (i - 1) * h;
        double k1 = h * f(y[i - 1], x);
        double k2 = h * f(y[i - 1] + 0.5 * k1, x + 0.5 * h);
        y[i] = y[i - 1] + k2;
    }
    return y;
}

// Компактная схема четвертого порядка
vector<double> fourthOrderScheme(double a, double b, double h, double y0) {
    int n = (b - a) / h + 1;
    vector<double> y(n);
    y[0] = y0;
    for (int i = 1; i < n; ++i) {
        double x = a + (i - 1) * h;
        double k1 = h * f(y[i - 1], x);
        double k2 = h * f(y[i - 1] + 0.5 * k1, x + 0.5 * h);
        double k3 = h * f(y[i - 1] + 0.5 * k2, x + 0.5 * h);
        double k4 = h * f(y[i - 1] + k3, x + h);
        y[i] = y[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    }
    return y;
}

// Расчет ошибок и локальных порядков сходимости
void calculateErrorsAndOrder(const vector<double>& y, double a, double h, ofstream& file) {
    int n = y.size();
    file << "j\txj\t\ty_exact\t\ty_num\t\tdelta\t\tp\n";
    file << fixed << setprecision(6);

    for (int i = 0; i < n; ++i) {
        double x = a + i * h;
        double y_exact = exactSolution(x);
        double delta = fabs(y_exact - y[i]);
        file << i << "\t" << x << "\t" << y_exact << "\t" << y[i] << "\t" << delta << "\t";
        if (i > 1) {
            // Локальный порядок сходимости
            double prevDelta = fabs(exactSolution(a + (i - 1) * h) - y[i - 1]);
            file << log(delta / prevDelta) / log(2.0) << "\n";
        } else {
            file << "N/A\n";
        }
    }
}

// Основная функция
int main() {
    double a, b, h;
    cout << "Введите значения a, b, h: ";
    cin >> a >> b >> h;

    double y0 = exactSolution(a); // Начальное условие

    // Разностные схемы
    auto y1 = firstOrderScheme(a, b, h, y0);
    auto y2 = secondOrderScheme(a, b, h, y0);
    auto y4 = fourthOrderScheme(a, b, h, y0);

    // Запись результатов в файл
    ofstream file("results.txt");

    file << "Разностная схема первого порядка:\n";
    calculateErrorsAndOrder(y1, a, h, file);

    file << "\nРазностная схема второго порядка:\n";
    calculateErrorsAndOrder(y2, a, h, file);

    file << "\nРазностная схема четвертого порядка:\n";
    calculateErrorsAndOrder(y4, a, h, file);

    file.close();
    cout << "Результаты сохранены в файл results.txt\n";

    return 0;
}
