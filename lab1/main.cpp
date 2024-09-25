#include "main.h"

#include <complex>
#include <iostream>
#include <optional>
#include <limits>
#include <string>
#include <vector>

double epsilon;

class quadric_equation {
private:
    const double a;
    const double b;
    const double c;
    double discriminant;
public:
    explicit quadric_equation(const double a, const double b, const double c) : a(a), b(b), c(c) {
        discriminant = b*b - 4*a*c;
    }

    std::vector<double> get_roots() const {
        std::vector<double> roots(2);

        if(discriminant < -epsilon) {

        } else if(discriminant > epsilon) {
            const double z1 = (-b + std::sqrt(discriminant)) / (2*a);
            const double z2 = (-b - std::sqrt(discriminant)) / (2*a);

            roots.push_back(z1);
            roots.push_back(z2);
        } else {
            const double z = (-b) / (2*a);
            roots.push_back(z);
        }

        return roots;
    }

};


class cubic_equation {
private:
    const double first;
    const double a;
    const double b;
    const double c;
public:
    cubic_equation(const double a, const double b, const double c, const double first = 1) : first(first), a(a), b(b), c(c){ }

    double calc_regular(const double x) const {
        return first*(x*x*x) + a*(x*x) + b*x + c;
    }

    double calc_derivative(const double x) const {
        return first*(3*x*x) + (2*a*x) + b;
    }

    quadric_equation get_derivative() const {
        return quadric_equation(first * 3, a * 2, b);
    }



    double get_a() const{
        return a;
    }

    double get_b() const{
        return b;
    }

    double get_c() const{
        return c;
    }

};

class bisection {
private:
    static const int step = 100;
public:
    static std::pair<double, double> find_interval(std::pair<double, double>& interval, const cubic_equation& f) {

        if(interval.first == std::numeric_limits<double>::infinity())

        const double c = 0;
        if((interval.first + interval.second) / 2)
    }
};

int main(int argc, char** argv) {
    if(argc != 5){
        std::cout << "Usage: " << argv[0] << " <epsilon> <a> <b> <c>" << std::endl;
        return 1;
    }

    double a, b, c;

    try{
        epsilon = std::stoi(std::string (argv[1]));
        a = std::stoi(std::string (argv[2]));
        b = std::stoi(std::string (argv[3]));
        c = std::stoi(std::string (argv[4]));
    } catch (const std::exception& e){
        std::cout << "Wrong argument!!! : " << e.what() << std::endl;
        return 1;
    }

    cubic_equation f(a, b, c);
    quadric_equation f_der = f.get_derivative();


    return 0;
}
