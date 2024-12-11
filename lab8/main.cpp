#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

typedef struct {
  int j;
  double x;
  double y1;
  double y2;
  double y3;
  double y_ex;
  double delta;
  double p;
} data_j;

double y(const double x) {
  return std::exp(std::sin(x));
}

double runge(const double yh1, const double yh2, const double yh3) {
  return std::log(std::abs(yh1 - yh2)/std::abs(yh2 - yh3)) / std::log(3);
}

std::vector<data_j> first_order(const double a, const double b, double h) {
  int n1 = static_cast<int>((b - a) / h) + 1;
  int n2 = static_cast<int>((b - a) / (h/3)) + 1;
  int n3 = static_cast<int>((b - a) / (h/9)) + 1;

  std::vector<data_j> result(n1);

  result[0] = {0, a, y(a), y(a), y(a), y(a), 0, 0};

  double y2[n2], y3[n3];

  y2[0] = y3[0] = y(a);

  for (int j = 1; j < n1; j++) {
    double x = a + j * h;
    result[j].j = j;
    result[j].x = x;
    result[j].y1 = result[j-1].y1 * (1 + h * std::cos(result[j-1].x));
  }
  for (int j = 1; j < n2; j++) {
    y2[j] = y2[j-1] * (1 + (h/3) * std::cos(a + (j-1) * (h/3)));
    if(j % 3 == 0) result[j/3].y2 = y2[j];
  }
  for (int j = 1; j < n3; j++) {
    y3[j] = y3[j-1] * (1 + (h/9) * std::cos(a + (j-1) * (h/9)));
    if(j % 9 == 0) result[j/9].y3 = y3[j];
  }

  for (int j = 1; j < n1; ++j) {
    result[j].y_ex = y(result[j].x);
    result[j].delta = std::abs(result[j].y_ex - result[j].y1);
    result[j].p = runge(result[j].y1, result[j].y2, result[j].y3);
  }

  return result;
}

std::vector<data_j> second_order(const double a, const double b, const double h) {
  int n1 = static_cast<int>((b - a) / h) + 1;
  int n2 = static_cast<int>((b - a) / (h/3)) + 1;
  int n3 = static_cast<int>((b - a) / (h/9)) + 1;

  std::vector<data_j> result(n1);

  result[0] = {0, a, y(a), y(a), y(a), y(a), 0, 0};

  double y2[n2], y3[n3];

  y2[0] = y3[0] = y(a);

  for (int j = 1; j < n1; j++) {
    double x = a + j * h;
    result[j].j = j;
    result[j].x = x;
    result[j].y1 = result[j-1].y1 * ((2 + h * cos(result[j-1].x))/(2 - h * std::cos(x)));
  }

  for (int j = 1; j < n2; j++) {
    double x = a + j * (h/3);
    y2[j] = y2[j-1] * ((2 + (h/3) * cos(a + (j-1) * (h/3)))/(2 - (h/3) * std::cos(x)));
    if(j % 3 == 0) result[j/3].y2 = y2[j];
  }

  for (int j = 1; j < n3; j++) {
    double x = a + j * (h/9);
    y3[j] = y3[j-1] * ((2 + (h/9) * cos(a + (j-1) * (h/9)))/(2 - (h/9) * std::cos(x)));
    if(j % 9 == 0) result[j/9].y3 = y3[j];
  }

  for (int j = 1; j < n1; ++j) {
    result[j].y_ex = y(result[j].x);
    result[j].delta = std::abs(result[j].y_ex - result[j].y1);
    result[j].p = runge(result[j].y1, result[j].y2, result[j].y3);
  }

  return result;
}

data_j fo(const data_j& zero, const double h) {
  data_j result;
  result.j = 1;
  result.x = zero.x + h;
  result.y1 = zero.y1 * (1 + h * std::cos(zero.x));
  result.y2 = zero.y2 * (1 + (h/3) * std::cos(zero.x));
  result.y3 = zero.y3 * (1 + (h/9) * std::cos(zero.x));
  result.y_ex = y(result.x);
  result.delta = std::abs(result.y_ex - result.y1);
  result.p = runge(result.y1, result.y2, result.y3);
  return result;
}

data_j so(const data_j& zero, const double h) {
  data_j result;
  result.j = 1;
  result.x = zero.x + h;
  result.y1 = zero.y1 * ((2 + h * cos(zero.x))/(2 - h * std::cos(result.x)));
  result.y2 = zero.y2 * ((2 + (h/3) * cos(zero.x))/(2 - (h/3) * std::cos(zero.x + (h/3))));
  result.y3 = zero.y3 * ((2 + (h/9) * cos(zero.x))/(2 - (h/9) * std::cos(zero.x + (h/9))));
  result.y_ex = y(result.x);
  result.delta = std::abs(result.y_ex - result.y1);
  result.p = runge(result.y1, result.y2, result.y3);
  return result;
}

data_j to(const data_j& zero, double h) {
  data_j result;
  result.j = 1;
  result.x = zero.x + h;
  result.y1 = zero.y1 * ((12 + 5*h*std::cos(zero.x)) / (h*std::cos(zero.x + 2*h)) - ((3 + h*std::cos(zero.x)) / (3 - h * std::cos(zero.x + 2*h))))
                / (((4*h*std::cos(zero.x + h)) / (3 - h*std::cos(zero.x + h*2))) + ((12 - 8*h*std::cos(zero.x + h)) / (h*std::cos(zero.x+h*2))));

  h /= 3;

  result.y2 = zero.y2 * ((12 + 5*h*std::cos(zero.x)) / (h*std::cos(zero.x + 2*h)) - ((3 + h*std::cos(zero.x)) / (3 - h * std::cos(zero.x + 2*h))))
                / (((4*h*std::cos(zero.x + h)) / (3 - h*std::cos(zero.x + h*2))) + ((12 - 8*h*std::cos(zero.x + h)) / (h*std::cos(zero.x+h*2))));

  h /= 3;

  result.y3 = zero.y3 * ((12 + 5*h*std::cos(zero.x)) / (h*std::cos(zero.x + 2*h)) - ((3 + h*std::cos(zero.x)) / (3 - h * std::cos(zero.x + 2*h))))
                / (((4*h*std::cos(zero.x + h)) / (3 - h*std::cos(zero.x + h*2))) + ((12 - 8*h*std::cos(zero.x + h)) / (h*std::cos(zero.x+h*2))));

  result.y_ex = y(result.x);
  result.delta = std::abs(result.y_ex - result.y1);
  result.p = runge(result.y1, result.y2, result.y3);
  return result;
}

std::vector<data_j> fourth_order(const double a, const double b, double h) {
  int n1 = static_cast<int>((b - a) / h) + 1;
  int n2 = static_cast<int>((b - a) / (h/3)) + 1;
  int n3 = static_cast<int>((b - a) / (h/9)) + 1;

  std::vector<data_j> result(n1);

  result[0] = {0, a, y(a), y(a), y(a), y(a), 0, 0};
  result[1] = to(result[0], h);

  double y2[n2], y3[n3];

  y2[0] = y3[0] = y(a);
  y2[1] = result[1].y2;
  y3[1] = result[1].y3;

  for (int j = 2; j < n1; j++) {
    double x = a + j * h;
    result[j].j = j;
    result[j].x = x;
    result[j].y1 = (result[j-2].y1*(3 + h * std::cos(a + (j-2) * h)) + 4 * h * result[j-1].y1 * std::cos(a + (j-1) * h))
                  /(3 - h * std::cos(x));
  }

  h /= 3;

  for (int j = 2; j < n2; j++) {
    double x = a + j * h;
    y2[j] = (y2[j-2]*(3 + h * std::cos(a + (j-2) * h)) + 4 * h * y2[j-1] * std::cos(a + (j-1) * h))
                  /(3 - h * std::cos(x));
    if(j % 3 == 0) result[j/3].y2 = y2[j];
  }

  h /= 3;

  for (int j = 2; j < n3; j++) {
    double x = a + j * h;
    y3[j] = (y3[j-2]*(3 + h * std::cos(a + (j-2) * h)) + 4 * h * y3[j-1] * std::cos(a + (j-1) * h))
                  /(3 - h * std::cos(x));
    if(j % 9 == 0) result[j/9].y3 = y3[j];
  }

  for (int j = 1; j < n1; ++j) {
    result[j].y_ex = y(result[j].x);
    result[j].delta = std::abs(result[j].y_ex - result[j].y1);
    result[j].p = runge(result[j].y1, result[j].y2, result[j].y3);
  }

  return result;
}

int main() {
  std::cout << "Input a, b, h: " << std::endl;
  double a, b, h;
  std::cin >> a >> b >> h;

  if (h <= 0 || a >= b) {
    std::cerr << "Invalid input" << std::endl;
    return -1;
  }

  std::ofstream out1("out1.csv");
  if(!out1) {
    std::cerr << "Error opening output file out.csv" << std::endl;
    return -1;
  }

  out1.clear();

  out1 << "j,x_j,y1_j,y2_j,y3_j,y_ex(x_j),delta_j,p_j\n";
  std::vector<data_j> fresult = first_order(a, b, h);

  for (auto & j : fresult) {
    out1 << j.j << "," << j.x << "," << j.y1 << "," << j.y2 << "," << j.y3 << "," << j.y_ex << "," << j.delta << "," << j.p <<std::endl;
  }

  std::ofstream out2("out2.csv");
  if(!out2) {
    std::cerr << "Error opening output file out.csv" << std::endl;
    return -1;
  }

  out2.clear();

  out2 << "j,x_j,y1_j,y2_j,y3_j,y_ex(x_j),delta_j,p_j\n";
  std::vector<data_j> sresult = second_order(a, b, h);

  for (auto & j : sresult) {
    out2 << j.j << "," << j.x << "," << j.y1 << "," << j.y2 << "," << j.y3 << "," << j.y_ex << "," << j.delta << "," << j.p <<std::endl;
  }

  std::ofstream out3("out3.csv");
  if(!out3) {
    std::cerr << "Error opening output file out.csv" << std::endl;
    return -1;
  }

  out3.clear();

  out3 << "j,x_j,y1_j,y2_j,y3_j,y_ex(x_j),delta_j,p_j\n";
  std::vector<data_j> tresult = fourth_order(a, b, h);

  for (auto & j : tresult) {
    out3 << j.j << "," << j.x << "," << j.y1 << "," << j.y2 << "," << j.y3 << "," << j.y_ex << "," << j.delta << "," << j.p <<std::endl;
  }

  return 0;
}
