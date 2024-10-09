#include <cstring>
#include <iostream>
#include <iomanip>

void fill_system1(long double *matrix, long double *f, const int N) {
  if (matrix == nullptr) {
    std::cout << "matrix is null" << std::endl;
    return;
  }

  std::memset(matrix, 0, N * N * sizeof(long double));
  std::memset(f, 0, N * sizeof(long double));

  for (int i = 0; i < N; ++i) {
    f[i] = 2;
    for (int j = 0; j < N; ++j) {
      if (i == j)
        matrix[i * N + j] = 2;
      if (std::abs(i - j) == 1)
        matrix[i * N + j] = -1;
    }
  }
}

void fill_system2(long double *matrix, long double *f, const int N,
                  const long double epsilon) {
  if (matrix == nullptr) {
    std::cout << "matrix is null" << std::endl;
    return;
  }

  std::memset(matrix, 0, N * N * sizeof(long double));
  std::memset(f, 0, N * sizeof(long double));

  for (int i = 0; i < N; ++i) {
    f[i] = 2 + epsilon;
    for (int j = 0; j < N; ++j) {
      if (i == j)
        matrix[i * N + j] = 2;
      if (std::abs(i - j) == 1)
        matrix[i * N + j] = -1;
    }
  }
}

void fill_system3(long double *matrix, long double *f, const int N, const long double gamma) {
  if (matrix == nullptr) {
    std::cout << "matrix is null" << std::endl;
    return;
  }

  std::memset(matrix, 0, N * N * sizeof(long double));
  std::memset(f, 0, N * sizeof(long double));

  for (int i = 0; i < N; ++i) {
    f[i] = 2 * (i + 2) + gamma;
    for (int j = 0; j < N; ++j) {
      if (i == j)
        matrix[i * N + j] = 2 * (i+1) + gamma;
      if (std::abs(i - j) == 1)
        matrix[i * N + j] = -1;
    }
  }
}

void print_matrix(const long double* matrix, const int N) {
  if (matrix == nullptr) {
    std::cout << "matrix is null" << std::endl;
    return;
  }

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      std::cout << matrix[i * N + j] << " ";
    }
    std::cout << std::endl;
  }
}

void print_roots(const long double* roots, const int N) {
  for (int i = 0; i < N; ++i) {
    std::cout << std::setprecision(9) << roots[i] << " ";
  }
}

void progonka(const long double *matrix, const long double *f,
                             long double* roots, const int N) {

  long double a_i, b_i, c_i, alpha[N], beta[N];

  b_i = -matrix[1];
  c_i = matrix[0];

  alpha[0] = b_i/c_i;
  beta[0] = f[0]/c_i;

  for (int i = 1; i < N; ++i) {
    a_i = -matrix[i*N + (i - 1)];
    b_i = -matrix[i*N + (i + 1)];
    c_i = matrix[i*N + i];

    alpha[i] = b_i/(c_i - a_i * alpha[i-1]);
    beta[i] = (f[i] + a_i * beta[i-1])/(c_i - a_i * alpha[i-1]);
  }

  roots[N - 1] = beta[N - 1];

  for (int i = N - 2; i >= 0; --i) {
    roots[i] = alpha[i] * roots[i + 1] + beta[i];
  }

}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cout << "Usage : " << argv[0] << " <N> <epsilon> <gamma>" << std::endl;
    return -1;
  }

  int N;
  long double epsilon, gamma;

  try {
    N = std::stoi(argv[1]);
    epsilon = std::stod(argv[2]);
    gamma = std::stod(argv[3]);
  } catch (const std::exception &e) {
    std::cout << "Wrong argument!!! : " << e.what() << std::endl;
    return -11;
  }

  auto *matrix = reinterpret_cast<long double *>(malloc(sizeof(long double) * N * N));
  auto *f = reinterpret_cast<long double *>(malloc(sizeof(long double) * N));
  long double roots[N];

  std::cout << "test 1" << std::endl;
  fill_system1(matrix, f, N);
  progonka(matrix, f, roots, N);
  print_roots(roots, N);
  std::cout << std::endl;

  std::cout << "test 2" << std::endl;
  fill_system2(matrix, f, N, epsilon);
  progonka(matrix, f, roots, N);
  print_roots(roots, N);
  std::cout << std::endl;

  std::cout << "test 3" << std::endl;
  fill_system3(matrix, f, N, gamma);
  progonka(matrix, f, roots, N);
  print_roots(roots, N);
  std::cout << std::endl;

  free(matrix);
  free(f);
  return 0;
}
