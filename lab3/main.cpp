#include <iostream>
#include <cstring>
#include <vector>

void fill_system1(double* matrix, double* f, const int N) {
    if(matrix == nullptr) {
        std::cout << "matrix is null" << std::endl;
        return;
    }

    std::memset(matrix, 0, N * N * sizeof(double));
    std::memset(f, 0, N * sizeof(double));

    for(int i = 0; i < N; ++i) {
        f[i] = 2;
        for(int j = 0; j < N; ++j) {
            if(i == j) matrix[i*N + j] = 2;
            if(std::abs(i - j) == 1) matrix[i*N + j] = -1;
        }
    }
}

void fill_system2(double* matrix, double* f, const int N, const double epsilon) {
    if(matrix == nullptr) {
        std::cout << "matrix is null" << std::endl;
        return;
    }

    std::memset(matrix, 0, N * N * sizeof(double));
    std::memset(f, 0, N * sizeof(double));

    for(int i = 0; i < N; ++i) {
        f[i] = 2 + epsilon;
        for(int j = 0; j < N; ++j) {
            if(i == j) matrix[i*N + j] = 2;
            if(std::abs(i - j) == 1) matrix[i*N + j] = -1;
        }
    }
}

void fill_system3(double* matrix, double* f, const int N, const double gamma) {
    if(matrix == nullptr) {
        std::cout << "matrix is null" << std::endl;
        return;
    }

    std::memset(matrix, 0, N * N * sizeof(double));
    std::memset(f, 0, N * sizeof(double));

    for(int i = 0; i < N; ++i) {
        f[i] = 2 * (i + 1) + gamma;
        for(int j = 0; j < N; ++j) {
            if(i == j) matrix[i*N + j] = 2 * i + gamma;
            if(std::abs(i - j) == 1) matrix[i*N + j] = -1;
        }
    }
}

void print_matrix(double* matrix, const int N) {
    if(matrix == nullptr) {
        std::cout << "matrix is null" << std::endl;
        return;
    }

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            std::cout << matrix[i*N + j] << " ";
        }
        std::cout << std::endl;
    }
}


std::vector<double> progonka(double* matrix, const int N) {
    std::vector<double> roots;

    double alpha_i, beta_i, a_i, b_i, c_i;

    b_i = 2;


    return roots;
}


int main(int argc, char** argv) {
    if(argc != 4) {
        std::cout << "Usage : " << argv[0] << " <N> <epsilon> <gamma>" << std::endl;
    }

    int N;
    double epsilon, gamma;

    try {
        N = std::stoi(argv[1]);
        epsilon = std::stod(argv[2]);
        gamma = std::stod(argv[3]);
    } catch (const std::exception& e) {
        std::cout << "Wrong argument!!! : " << e.what() << std::endl;
        return 1;
    }

    auto* matrix = reinterpret_cast<double*>(malloc(sizeof(double) * N * N));
    auto* f = reinterpret_cast<double*>(malloc(sizeof(double) * N));


    print_matrix(matrix, N);

    free(matrix);
    free(f);
    return 0;
}
