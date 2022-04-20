#include "extra.h"


SLAR::SLAR()
{
  order = 3;
  AMatrix = new double*[3];
  for(int i{0}; i < 3; i++)
  {
    AMatrix[i] = new double[3];
  }
  BMatrix = new double[3];
  XMatrix = new double [3];
}

SLAR::SLAR(int orderM)
  : order{orderM}
{
  AMatrix = new double*[orderM];
  for(int i{0}; i < orderM; i++)
  {
    AMatrix[i] = new double[orderM];
  }
  BMatrix = new double[orderM];
  XMatrix = new double[orderM];
}

void SLAR::print()
{
  std::string divideLine(80, '-');
  std::cout << divideLine << std::endl;
  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      if(AMatrix[i][j] >= 0)
      {
        std::cout << " + " << AMatrix[i][j] << "x[" << j + 1 << "] ";
      } else
      {
        std::cout << " - " << AMatrix[i][j] * (-1) << "x[" << j + 1 << "] ";
      }
    }
    std::cout << " = " << BMatrix[i] << std::endl;
  }
  std::cout << divideLine << std::endl;
}

void SLAR::printExtendedMatrix()
{
  std::string divideLine(80, '-');
  std::cout << divideLine << std::endl;
  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      std::cout << AMatrix[i][j] << " ";
    }
    std::cout << "| " << BMatrix[i] << std::endl;
  }
  std::cout << divideLine << std::endl;
}

std::istream& operator>>(std::istream& in, SLAR system)
{
  for(int i{0}; i < system.order; i++)
  {
    for(int j{0}; j < system.order; j++)
    {
      std::cout << "Enter element of Matrix A in " << i+1 << " row and " << j+1 << " column" << std::endl;
      in >> system.AMatrix[i][j];
    }
  }

  for(int i{0}; i < system.order; i++)
  {
    std::cout << "Enter element of Matrix B in " << i+1 << " row" << std::endl;
    in >> system.BMatrix[i];
  }
}

void SLAR::readDataFromFile(std::string filepath)
{
  std::ifstream file(filepath);

  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      file >> AMatrix[i][j];
    }
    file >> BMatrix[i];
  }
  file.close();
}

void SLAR::calcMinor(double **matrix, double **minor, int col, int row, int orderM)
{
  int ki, kj, di, dj;
  di = 0;
  for (ki = 0; ki < orderM - 1; ki++)
  {
    if (ki == col) di = 1;
    dj = 0;
    for (kj = 0; kj < orderM - 1; kj++)
    {
      if (kj == row) dj = 1;
      minor[ki][kj] = matrix[ki + di][kj + dj];
    }
  }
}

double SLAR::calcDeterminant(double **matrix, int orderM)
{
  int i, j, k, n;
  double det;
  double **p;
  p = new double*[orderM];
  for (i = 0; i < orderM; i++)
    p[i] = new double[orderM];
  j = 0; det = 0;
  k = 1;
  n = orderM - 1;
  if (orderM<1) std::cout << "Impossible to calc determinant" << std::endl;
  if (orderM == 1) {
    det = matrix[0][0];
    return(det);
  }
  if (orderM == 2) {
    det = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
    return(det);
  }
  if (orderM>2) {
    for (i = 0; i < orderM; i++) {
      calcMinor(matrix, p, i, 0, orderM);
      det = det + k * matrix[i][0] * calcDeterminant(p, n);
      k = -k;
    }
  }
  return(det);
}

void SLAR::swapRows(int row1, int row2)
{
  double temp[order+1];

  for(int i{0}; i < order; i++)
  {
    temp[i] = AMatrix[row1][i];
  }
  temp[order] = BMatrix[row1];
  for(int i{0}; i < order; i++)
  {
    AMatrix[row1][i] = AMatrix[row2][i];
  }
  BMatrix[row1] = BMatrix[row2];
  for(int i{0}; i < order; i++)
  {
    AMatrix[row2][i] = temp[i];
  }
  BMatrix[row2] = temp[order];
}

void SLAR::swapCols(int col1, int col2)
{
  double temp[order+1];

  for(int i{0}; i < order; i++)
  {
    temp[i] = AMatrix[i][col1];
  }
  temp[order] = BMatrix[col1];
  for(int i{0}; i < order; i++)
  {
    AMatrix[i][col1] = AMatrix[i][col2];
  }
  for(int i{0}; i < order; i++)
  {
    AMatrix[i][col2] = temp[i];
  }
}

void SLAR::addRow(int rowToAdd, int rowAdded, double mult)
{
  for(int i{0}; i < order; i++)
  {
    AMatrix[rowToAdd][i] += AMatrix[rowAdded][i] * mult;
  }
  BMatrix[rowToAdd] += BMatrix[rowAdded] * mult;
}

void SLAR::GaussianMethod()
{
  std::cout << "************************* Gaussian Method **************************************"
    << std::endl;
  printExtendedMatrix();
  double biggest{0};
  int mainRow{0};
  for(int i{0}; i < order-1; i++)
  {
    biggest = AMatrix[i][i];
    mainRow = i;
    for(int j{i}; j < order; j++)
    {
      if(biggest < fabs(AMatrix[j][i]))
      {
        biggest = fabs(AMatrix[j][i]);
        mainRow = j;
      }
    }
    if(biggest == 0)
    {
      std::cout << "Impossible to solve this system" << std::endl;
      return;
    }
    if(mainRow != i)
    {
      swapRows(mainRow, i);
    }
    double mult{0};
    for(int j{i+1}; j < order; j++)
    {
      mult = AMatrix[j][i] / AMatrix[i][i];
      addRow(j, i, -mult);
    }
    printExtendedMatrix();
  }

  XMatrix[order-1] = BMatrix[order-1]/AMatrix[order-1][order-1];

  for(int i{1}; i < order; i++)
  {
    double sub{0};
    for(int j{0}; j < order; j++)
    {
      sub+= XMatrix[order-1-j]*AMatrix[order-1-i][order-1-j];
    }
    XMatrix[order-1-i] = (BMatrix[order-1-i]-sub)/AMatrix[order-1-i][order-1-i];
  }


  std::cout << std::endl;
  for(int i{0}; i < order; i++)
  {
    std::cout << "X[" << i+1 << "] = " << XMatrix[i] << std::endl;
  }
}

void SLAR::LU_DecompositionMethod()
{
  std::cout << "************************* LU-Decomposition Method ******************************"
            << std::endl;

  double matrL[order][order], matrU[order][order];

  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      matrL[i][j] = 0;
      matrU[i][j] = 0;
    }
  }

  for(int i{0}; i < order; i++)
  {
    matrL[i][0] = AMatrix[i][0];
    matrU[i][i] = 1;
    matrU[0][i] = AMatrix[0][i]/AMatrix[0][0];
  }

  for(int i{1}; i < order; i++)
  {
    for(int j{1}; j < i+1; j++)
    {
      double sub{0};
      for(int k{0}; k < i; k++)
      {
        sub+= matrL[i][k]*matrU[k][j];
      }
      matrL[i][j] = AMatrix[i][j] - sub;
    }
    for(int j{i+1}; j < order; j++)
    {
      double sub{0};
      for(int k{0}; k < i; k++)
      {
        sub+= matrL[i][k]*matrU[k][j];
      }
      matrU[i][j] = (AMatrix[i][j] - sub)/matrL[i][i];
    }
  }

  std::string divideLine(80, '-');
  std::cout << divideLine << std::endl;
  std::cout << "Matrix L" << std::endl;
  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      std::cout << matrL[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << divideLine << std::endl;
  std::cout << divideLine << std::endl;
  std::cout << "Matrix U" << std::endl;
  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      std::cout << matrU[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << divideLine << std::endl;

  double vectorY[order];

  vectorY[0] = BMatrix[0]/matrL[0][0];

  for(int i{1}; i < order; i++)
  {
    double sub{0};
    for(int j{0}; j < i; j++)
    {
      sub += matrL[i][j] * vectorY[j];
    }
    vectorY[i] = (BMatrix[i] - sub)/matrL[i][i];
  }

  std::cout << divideLine << std::endl;
  for(int i{0}; i < order; i++)
  {
    std::cout << "Y[" << i+1 << "] = " << vectorY[i] << " ";
  }
  std::cout << std::endl;
  std::cout << divideLine << std::endl;

  XMatrix[order-1] = vectorY[order-1];

  for(int i{1}; i < order; i++)
  {
    double sub{0};
    for(int j{0}; j < order; j++)
    {
      sub += matrU[order-1-i][order-1-j] * XMatrix[order-1-j];
    }
    XMatrix[order-1-i] = vectorY[order-1-i] - sub;
  }
  std::cout << divideLine << std::endl;
  for(int i{0}; i < order; i++)
  {
    std::cout << "X[" << i+1 << "] = " << XMatrix[i] << " ";
  }
  std::cout << std::endl;
  std::cout << divideLine << std::endl;
}