#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

typedef struct {
  int j;
  double x_j;
  double y_xj;
  double y_ex_xj;
  double delta_j;
  double p_j;
} data;

double y(const double x) {
  return std::exp(std::sin(x));
}

std::vector<data> first_order(const double a, const double b, const double h) {
  std::vector<data> result;
  const int n = static_cast<int>((b - a) / h);
  result.reserve(n);

  data cur_data = {0, a, y(a), y(a), 0, 0};
  result.push_back(cur_data);

  for (int j = 1; j < n; ++j) {
    cur_data.j = j;
    cur_data.x_j = a + h * j;
    cur_data.y_ex_xj = result[j-1].y_xj + h * result[j-1].x_j * std::cos(cur_data.x_j);
    cur_data.delta_j = cur_data.y_ex_xj - cur_data.y_xj;
    result.push_back(cur_data);
  }
  return result;
}

int main() {
  std::cout << "Input a, b, h: " << std::endl;
  double a, b, h;
  std::cin >> a >> b >> h;

  std::cout << h << std::endl;

  if (h <= 0 || a >= b) {
    std::cerr << "Invalid input" << std::endl;
    return -1;
  }

  std::ofstream out("out.csv");
  if(!out) {
    std::cerr << "Error opening output file out.csv" << std::endl;
    return -1;
  }

  out << "j,x_j,y_ex(x_j),delta_j,p_j\n";
  std::vector<data> result = first_order(a, b, h);

  for (auto & j : result) {
    out << j.j << "," << j.x_j << "," << j.y_ex_xj << "," << j.delta_j << "," << j.p_j <<std::endl;
  }

  return 0;
}
