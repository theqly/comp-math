#include "bisection.h"

#include <iostream>

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cout << "Usage: " << argv[0] << " <epsilon> <a> <b> <c>" << std::endl;
    return 1;
  }

  double a, b, c;

  try {
    epsilon = std::stod(argv[1]);
    a = std::stod(argv[2]);
    b = std::stod(argv[3]);
    c = std::stod(argv[4]);
  } catch (const std::exception &e) {
    std::cout << "Wrong argument!!! : " << e.what() << std::endl;
    return 1;
  }

  const cubic_equation f(a, b, c);

  /*std::vector<double> roots = f.get_roots();

  // Вывод корней
  for (const auto &root : roots) {
    std::cout << "Root: " << root << std::endl;
  }*/

  std::cout << "root: " << bisection::find_root(f, -2, -1).value() << std::endl;


  return 0;
}
