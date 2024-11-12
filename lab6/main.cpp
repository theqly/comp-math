#include <iostream>
#include <iomanip>
#include <cmath>

const double a = 5.0;
const double b = 7.0;

double f(const double x) {
    return std::exp(x) * std::cos(x);
}

double f2(const double x) {
    return -2 * std::exp(x) * std::sin(x);
}

double f4(const double x) {
    return -4 * std::exp(x) * std::cos(x);
}

double runge(const double S1, const double S2, const double S3) {
    const double frac = (S1 - S2) / (S2 - S3);
    return std::log2(std::fabs(frac));
}

double exact_value() {
    return (std::exp(7) * (std::sin(7) + std::cos(7)) - std::exp(5) * (std::sin(5) + std::cos(5))) / 2;
}

double trapezoid_method(const int N, const double delta) {
    double S = 0.0;
    const double h = delta;
    double xi = a;
    for(int i = 0; i < N; ++i) {
        const double Si = ((f(xi) + f(xi+h))/2)*h;
        xi += h;
        S += Si;
    }
    return S;
}

double parabola_method(const int N, const double delta) {
    double S = 0.0;
    const double h = delta / 2;

    for(int i = 0; i < N; ++i) {
        double xi = a + i*delta;

        const double Si = h * (f(xi) + f(xi + 2*h) + 4*f(xi +h)) / 3;
        S += Si;
    }
    return S;
}

double quadrature_method(const int N, const double delta) {
    double S = 0.0;
    const double h = delta / 6;
    for(int i = 0; i < N; ++i) {
        double xi = a + i*delta;
        const double Si = ((3 * h) / 4) * (f(xi) + 3*(f(xi+2*h) + f(xi+4*h)) + f(xi+6*h));
        S += Si;
    }
    return S;
}


int main(int argc, char** argv) {
    if(argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return -1;
    }

    int N = std::stoi(argv[1]);
    double delta = (b - a)/N;

    std::cout << std::fixed << std::setprecision(15);

    double I = exact_value();
    std::cout << "точное значение: I = " << I << std::endl;

    double St1 = trapezoid_method(N, delta);
    std::cout << "метод трапеций: S = " << St1 << std::endl;

    double Sp1 = parabola_method(N, delta);
    std::cout << "метод парабол: S = " << Sp1 << std::endl;

    double Sq1 = quadrature_method(N, delta);
    std::cout << "метод квадратурной формулы: S = " << Sq1 << std::endl;

    //погрешности

    std::cout << "-----ПОГРЕШНОСТИ-----" << std::endl;

    std::cout << "метод трапеций: " << I - St1 << std::endl;

    std::cout << "метод парабол: " << I - Sp1 << std::endl;

    std::cout << "метод квадратурной формулы: " << I - Sq1 << std::endl;


    //рунге

    N *= 2;
    delta = (b - a)/N;

    double St2 = trapezoid_method(N, delta);

    double Sp2 = parabola_method(N, delta);

    double Sq2 = quadrature_method(N, delta);

    N *= 2;
    delta = (b - a)/N;

    double St3 = trapezoid_method(N, delta);

    double Sp3 = parabola_method(N, delta);

    double Sq3 = quadrature_method(N, delta);

    std::cout << "----- RUNGE -----" << std::endl;

    std::cout << "метод трапеций: S = " << runge(St1, St2, St3) << std::endl;

    std::cout << "метод парабол: S = " << runge(Sp1, Sp2, Sp3) << std::endl;

    std::cout << "метод квадратурной формулы: S = " << runge(Sq1, Sq2, Sq3) << std::endl;

    return 0;
}
