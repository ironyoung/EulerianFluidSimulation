#include <iostream>
#include "FLUID_3D.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  cout << "=========================================================================" << endl;
  cout << " Physically-based Eulerian Fluid Simulation: Smoke Simulator " << endl;
  cout << "=========================================================================" << endl;

  int xRes = 48;
  int yRes = 64;
  int zRes = 48;
  int amplify = 1;
  xRes *= amplify;
  yRes *= amplify;
  zRes *= amplify;

  int totalCells = xRes * yRes * zRes;
  int amplifiedCells = totalCells * amplify * amplify * amplify;
  
  // print out memory requirements
  long long int coarseSize = sizeof(double) * totalCells * 22 + 
                   sizeof(unsigned char) * totalCells;
  long long int totalMB = coarseSize / 1048576;
  cout << " Current coarse resolution: " << xRes << " x " << yRes << " x " << zRes << endl;
  cout << " At least " << totalMB << " MB of RAM needed " << endl;
  cout << "=========================================================================" << endl;
  cout.flush();

  // create output directories
  system("mkdir result");
  
  // 声明，即调用了构造函数
  FLUID_3D fluid(xRes, yRes, zRes, amplify);
  for (int x = 0; x < 100; x++)
  {
    fluid.addSmokeColumn();
    fluid.step();
  }

	return EXIT_SUCCESS;
}
