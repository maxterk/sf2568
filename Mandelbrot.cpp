#include <iostream>
#include <mpi.h>
#include <cmath>

using namespace std;

class Pixel
{
    double re,im;
    int pixelValue=-1;
public:
  Pixel(double a, double b)
  {
    re=a;
    im=b;
  }
  int calcPixel(double b, int N)
  {
    if (pixelValue==-1) {
      //This is z_0
      double tempReal=0, tempImaginary=0,
        nextReal,nextImaginary,currentAbsolute;
      int iters=1;
      double bb=pow(b,2);
      do
      {
        iters++;
        nextReal= re+pow(tempReal,2)-pow(tempImaginary, 2);
        nextImaginary=im+2*tempReal*tempImaginary;
        currentAbsolute=pow(nextReal,2)+pow(nextImaginary,2);
        tempReal=nextReal;
        tempImaginary=nextImaginary;
      } while(currentAbsolute<bb && iters<N);
      pixelValue=iters;
      return pixelValue;
    }
    else
    {
      return pixelValue;
    }

  }

};
int main(int argc, char *argv[]) {
  int rc = MPI_Init(&argc,&argv);
  int size,rank;
  MPI_Status status;

  int hello = MPI_Comm_size(MPI_COMM_WORLD, &size);
  int hello2 = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Pixel p(1/(1+(float) rank),1/(1+(float) rank));
  std::cout <<"Rank: "<< rank <<'\n'<<p.calcPixel(2.0,256)<< "\n\n";

  MPI_Finalize();
  return 0;


}
