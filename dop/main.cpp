#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

const double h = 0.1;
const double r = 0.8;
const double T = 4.0;
const double a = -1.0, b = 10.0;

const double tau = r * h;

double initial_condition(double x) { return (x < 0) ? 1.0 : 0.0; }

std::vector<std::vector<double>> exact_solution(const std::vector<double> &x, const std::vector<double> &t) {

  std::vector<std::vector<double>> u(t.size(),
                                     std::vector<double>(x.size(), 0.0));
  for(size_t n = 0; n < u.size(); ++n) {
    for (size_t i = 0; i < x.size(); ++i) {
      u[n][i] = initial_condition(x[i] - t[n]);
    }
  }
  return u;
}

std::vector<std::vector<double>> godunov_explicit(const std::vector<double> &x, const std::vector<double> &t) {

  std::vector<std::vector<double>> u(t.size(),
                                     std::vector<double>(x.size(), 0.0));
  for (int j = 0; j < x.size(); ++j) {
    u[0][j] = initial_condition(x[j]);
  }

  for (int n = 0; n < u.size() - 1; ++n) {
    u[n + 1][0] = initial_condition(a);

    for (int j = 1; j < x.size(); ++j) {
      u[n + 1][j] = u[n][j] * (1 - r) + r * u[n][j - 1];
    }

    u[n + 1][x.size() - 1] = initial_condition(b);
  }

  return u;
}

std::vector<std::vector<double>> implicit_scheme(const std::vector<double> &x, const std::vector<double> &t) {
  std::vector<std::vector<double>> u(t.size(), std::vector<double>(x.size(), 0.0));

  for (size_t i = 0; i < x.size(); ++i) {
    u[0][i] = initial_condition(x[i]);
  }

  std::vector<std::vector<double>> A(x.size(), std::vector<double>(x.size(), 0.0));
  for (size_t i = 1; i < x.size() - 1; ++i) {
    A[i][i - 1] = -r / 2.0;
    A[i][i] = 1.0;
    A[i][i + 1] = r / 2.0;
  }
  A[0][0] = 1.0;
  A[x.size() - 1][x.size() - 1] = 1.0;

  std::vector<double> B(x.size());
  for (size_t n = 0; n < t.size() - 1; ++n) {
    for (size_t i = 1; i < x.size() - 1; ++i) {
      B[i] = u[n][i];
    }
    B[0] = initial_condition(x[0] - t[n + 1]);
    B[x.size() - 1] = initial_condition(x.back() - t[n + 1]);

    std::vector<std::vector<double>> A_temp = A;
    std::vector<double> y = B;

    for (size_t i = 1; i < x.size(); ++i) {
      double factor = A_temp[i][i - 1] / A_temp[i - 1][i - 1];
      for (size_t j = i - 1; j <= i + 1 && j < x.size(); ++j) {
        A_temp[i][j] -= factor * A_temp[i - 1][j];
      }
      y[i] -= factor * y[i - 1];
    }

    u[n + 1][x.size() - 1] = y[x.size() - 1] / A_temp[x.size() - 1][x.size() - 1];
    for (int i = x.size() - 2; i >= 0; --i) {
      u[n + 1][i] = (y[i] - A_temp[i][i + 1] * u[n + 1][i + 1]) / A_temp[i][i];
    }
  }

  return u;
}

void plot_results(const std::vector<double> &x,
                  const std::vector<std::vector<double>> &u_exact,
                  const std::vector<std::vector<double>> &u_godunov,
                  const std::vector<std::vector<double>> &u_implicit,
                  const std::vector<double> &t) {
  std::ofstream file("results.dat");


  for (size_t t_idx = 0; t_idx < t.size(); ++t_idx) {
    file << "# Time = " << t[t_idx] << std::endl;
    file << "# x Exact Godunov Implicit" << std::endl;
    for (size_t i = 0; i < x.size(); ++i) {
      file << x[i] << " " << u_exact[t_idx][i] << " " << u_godunov[t_idx][i]
           << " " << u_implicit[t_idx][i] << std::endl;
    }
    file << std::endl;
  }

  file.close();
}

int main() {

  int x_steps = static_cast<int>((b - a) / h) + 1;
  int t_steps = static_cast<int>(T / tau) + 1;

  std::vector<double> x;
  for (int j = 0; j < x_steps; ++j) {
    double xj = a + j * h;
    x.push_back(xj);
  }

  std::vector<double> t;
  for (int n = 0; n < t_steps; ++n) {
    double tn = n * tau;
    t.push_back(tn);
  }

  auto u_exact = exact_solution(x, t);
  auto u_godunov = godunov_explicit(x, t);
  auto u_implicit = implicit_scheme(x, t);

  plot_results(x, u_exact, u_godunov, u_implicit, t);

  return 0;
}
