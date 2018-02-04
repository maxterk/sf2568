#include <iostream>
#include <mpi.h>
#include <cmath>

using namespace std;

class Pixel
{
    double re,im;
    int pixelValue=-1;
public:
  Pixel(double a=0, double b=0)
  {
    re=a;
    im=b;
  }
  
  Pixel set_pixel(double a, double b){
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

  double b=2;
  int w=2048;
  int h=2048;
  int r=w*h;
  double dx= (double) 2*b/(w-1);
  double dy= (double) 2*b/(h-1);
  double dreal;
  double dimag;
  FILE *fp;
  fp = fopen("color.txt","w");
  int rc = MPI_Init(&argc,&argv);
  int size,rank;
  MPI_Status status;
  int hello = MPI_Comm_size(MPI_COMM_WORLD, &size);
  int hello2 = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  
  
 Pixel * pixel_array;
 unsigned char* color;
  
  pixel_array= new Pixel [w*h];
  color=new unsigned char [w*h];

  for(int x=0;x<w;x++){
    dreal=(double)x*dx-b;
    for(int y=0;y<h;y++ ){
      dimag=(double)y*dy-b;
      pixel_array[x*h+y].set_pixel(dreal,dimag);     
      //cout<< pixel_array[x*h+y].calcPixel(2.0,256)<< "\n";
    fprintf(fp, "%hhu", pixel_array[x*h+y].calcPixel(2.0,256));
    fprintf(fp, "\n");
    }

    }

 
  
  MPI_Finalize();
  fclose(fp);
  return 0;


}
