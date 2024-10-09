#include "bisection.h"

#include <complex>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

//-----------------------------------QUADRIC-----------------------------------

quadric_equation::quadric_equation(const double a, const double b,
                                   const double c)
    : a(a), b(b), c(c) {
  discriminant = b * b - 4 * a * c;
}

double quadric_equation::calc(const double x) const {
  return a * (x * x) + b * x + c;
}

std::vector<double> quadric_equation::get_roots() const {

  std::vector<double> roots;

  if (discriminant < -epsilon) {
  } else if (discriminant > epsilon) {
    const double z1 = (-b + std::sqrt(discriminant)) / (2 * a);
    const double z2 = (-b - std::sqrt(discriminant)) / (2 * a);
    if (z1 < z2) {
      roots.push_back(z1);
      roots.push_back(z2);
    } else {
      roots.push_back(z2);
      roots.push_back(z1);
    }

  } else {
    const double z = (-b) / (2 * a);
    roots.push_back(z);
  }

  return roots;
}

//-----------------------------------CUBIC-----------------------------------

 cubic_equation::cubic_equation(const double a, const double b, const double c,
                               const double first) : first(first), a(a), b(b), c(c) {}




double cubic_equation::calc(const double x) const {
  return first * (x * x * x) + a * (x * x) + b * x + c;
}

quadric_equation cubic_equation::get_derivative() const {
  return quadric_equation{first * 3, a * 2, b};
}

std::vector<double> cubic_equation::get_roots() const {
  std::vector<double> roots;

  quadric_equation derivative = get_derivative();
  std::vector<double> derivative_roots = derivative.get_roots();

  switch (derivative_roots.size()) {
  case 0: {

    std::cout << "0 roots in derivative" << std::endl;

    double result = calc(0);

    if (result < -epsilon) {
      const double root = bisection::find_root(
                              *this, 0, std::numeric_limits<double>::infinity())
                              .value_or(-666);
      roots.push_back(root);
    } else if (result > epsilon) {
      const double root = bisection::find_root(
                              *this, std::numeric_limits<double>::infinity(), 0)
                              .value_or(-666);
      roots.push_back(root);
    } else {
      roots.push_back(0);
    }

    break;

  }
  case 1: {

    std::cout << "1 root in derivative: " << derivative_roots[0] << std::endl;

    //const double zyuzya = calc(-(derivative_roots[0]/3));
    const double zyuzya = calc(-(a/3));
    double root;

    if (std::abs(zyuzya) <= epsilon) {
      root = -(a/3);
      roots.push_back(root);
      roots.push_back(root);
      roots.push_back(root);
    } else if(zyuzya < -epsilon){
      root = bisection::find_root(*this, zyuzya, std::numeric_limits<double>::infinity()).value();
      roots.push_back(root);
      roots.push_back(root);
      roots.push_back(root);
    } else {
      root = bisection::find_root(*this, std::numeric_limits<double>::infinity(), zyuzya).value();
      roots.push_back(root);
      roots.push_back(root);
      roots.push_back(root);
    }

    break;
  }
  case 2: {

    std::cout << "2 roots in derivative: " << derivative_roots[0] << ", " << derivative_roots[1] << std::endl;

    const double alpha = derivative_roots[0];
    const double beta = derivative_roots[1];

    const double f_alpha = calc(alpha);
    const double f_beta = calc(beta);

    // BIGGEST IF IN MY LIFE

    if (f_alpha > epsilon && f_beta < -epsilon) {

      std::cout << "case 1" << std::endl;

      const double root1 = bisection::find_root(*this, std::numeric_limits<double>::infinity(), alpha).value_or(-666);
      const double root2 = bisection::find_root(*this, alpha, beta).value_or(-666);
      const double root3 = bisection::find_root(*this, beta, std::numeric_limits<double>::infinity()).value_or(-666);
      roots.push_back(root1);
      roots.push_back(root2);
      roots.push_back(root3);
    } else if (f_alpha > epsilon && std::abs(f_beta) < epsilon) {

      std::cout << "case 2" << std::endl;

      const double root = bisection::find_root(*this, std::numeric_limits<double>::infinity(), alpha).value_or(-666);
      roots.push_back(beta);
      roots.push_back(beta);
      roots.push_back(root);
    } else if (f_beta < -epsilon && std::abs(f_alpha) < epsilon) {

      std::cout << "case 3" << std::endl;

      const double root = bisection::find_root(*this, beta, std::numeric_limits<double>::infinity()).value_or(-666);
      roots.push_back(alpha);
      roots.push_back(alpha);
      roots.push_back(root);
    } else if (f_alpha > epsilon && f_beta > epsilon) {

      std::cout << "case 4" << std::endl;

      const double root = bisection::find_root(*this, std::numeric_limits<double>::infinity(), alpha).value_or(-666);
      roots.push_back(root);
    } else if(f_alpha < -epsilon && f_beta < -epsilon) {

      std::cout << "case 5" << std::endl;

      const double root = bisection::find_root(*this, f_beta, std::numeric_limits<double>::infinity()).value_or(-666);
      roots.push_back(root);
    } else {
      const double root = (alpha + beta) / 2;
      roots.push_back(root);
      roots.push_back(root);
      roots.push_back(root);
    }

    break;
  }
  }

  return roots;
}

//-----------------------------------BISECTION-----------------------------------

std::optional<double> bisection::find_root(const ::cubic_equation &f,
                                           double left, double right) {
  if (left == std::numeric_limits<double>::infinity())
    left = MIN_TRY;
  if (right == std::numeric_limits<double>::infinity())
    right = MAX_TRY;

  if (f.calc(left) * f.calc(right) > 0) {
    std::cerr << "no root in interval: [" << left << ", " << right << "] "
              << std::endl;
    return std::nullopt;
  }

  for (long long int i = 0; i < MAX_ITER; ++i) {
    double mid = (left + right) / 2.0;
    double f_mid = f.calc(mid);

    if (std::abs(f_mid) < epsilon) {
      std::cout << "iterations: " << i << std::endl;
      return mid;
    }

    if (f.calc(left) * f_mid < 0) {
      right = mid;
    } else {
      left = mid;
    }
  }

  return (left + right) / 2.0;
}
