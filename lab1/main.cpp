#include <iostream>
#include <string>
#include <cmath>

double epsilon;

class function {
private:
    const double a;
    const double b;
    const double c;
public:
    function(const double a, const double b, const double c) : a(a), b(b), c(c) { }

    double regular(const double x) const {
        return (x*x*x) + a*(x*x) + (b*x) + c;
    }

    double derivative(const double x) const {
        return (3*x*x) + (2*a*x) + b;
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

class equation {
private:
    const function f;
    double discriminant;
public:
    explicit equation(const function& f) : f(f){
        discriminant = f.get_b()*f.get_b() - 4 * f.get_a() * f.get_c();
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
    }




    return 0;
}
