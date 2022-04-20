#include "extra.h"

int main()
{
  int order{3};

//  std::cout << "Enter matrix A order: " << std::endl;
//  std::cin >> order;
//  while(order < 1)
//  {
//    std::cout << "You entered wrong order, please try again" << std::endl;
//    std::cin >> order;
//  }
  SLAR MySystem1(order);
  SLAR MySystem2(order);

  MySystem1.readDataFromFile("D:\\Work\\University\\C1S2\\Numerical Methods\\lab04\\system.txt");
  MySystem2.readDataFromFile("D:\\Work\\University\\C1S2\\Numerical Methods\\lab04\\system.txt");
//  std::cin >> MySystem1;
//  std::cin >> MySystem2;
//  MySystem1.GaussianMethod();

  MySystem1.print();
  MySystem1.GaussianMethod();
  MySystem2.LU_DecompositionMethod();
  return 0;
}
