#ifndef LAB04_EXTRA_H
#define LAB04_EXTRA_H
#include <iostream>
#include <cmath>
#include <fstream>
class SLAR
{
public:
  SLAR();
  SLAR(int orderM);
private:
  double **AMatrix;
  double *BMatrix;
  double *XMatrix;
  int order;
public:
  double **getAMatrix() {return AMatrix;}
  double *getBMatrix() {return BMatrix;}
  int getOrder() {return order;}
  void readDataFromFile(std::string filepath);
  void print();
  void printExtendedMatrix();
  void calcMinor(double **matrix, double **minor, int row, int col, int orderM);
  double calcDeterminant(double **matrix, int orderM);
  void swapRows(int row1, int row2);
  void swapCols(int col1, int col2);
  void addRow(int rowToAdd, int rowAdded, double mult);
  void GaussianMethod();
  void LU_DecompositionMethod();

  friend std::istream& operator>>(std::istream& in, SLAR system);
};

std::istream& operator>>(std::istream& in, SLAR system);
#endif //LAB04_EXTRA_H
