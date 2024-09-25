#ifndef BISECTION_H
#define BISECTION_H
#include <optional>
#include <vector>

//NEED TO BE DEFINED FROM MAIN !!! !!! !!!
inline double epsilon;

class quadric_equation {
private:
  const double a;
  const double b;
  const double c;
  double discriminant;

public:
  quadric_equation(double a, double b, double c);

  double calc(double x) const;

  std::vector<double> get_roots() const;
};

class cubic_equation {
private:
  const double first;
  const double a;
  const double b;
  const double c;

public:
  cubic_equation(const double a, const double b, const double c,
                 const double first = 1);

  double calc(double x) const;

  quadric_equation get_derivative() const;

  std::vector<double> get_roots() const;
};

class bisection {
private:
  static const long long int MAX_ITER = 10000000000;

  static const int MIN_TRY = -1000000;
  static const int MAX_TRY = 1000000;

public:
  static std::optional<double> find_root(const cubic_equation &f, double left,
                                         double right);
};

#endif // BISECTION_H
