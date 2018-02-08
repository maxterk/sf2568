# include <iostream>
# include <mpi.h>
# include <cmath>

using namespace std;

int calcPix(double maxAbsoluteValue, int maxIter, double re, double im)
{
  double tempReal=0, tempImaginary=0,
  nextReal,nextImaginary,currentAbsolute;
  int iters=0;
  double maxAbsoluteValueSquared=pow(maxAbsoluteValue,2);
  do
  {
    iters++;
    nextReal= re+pow(tempReal,2)-pow(tempImaginary, 2);
    nextImaginary=im+2*tempReal*tempImaginary;
    currentAbsolute=pow(nextReal,2)+pow(nextImaginary,2);
    tempReal=nextReal;
    tempImaginary=nextImaginary;
  } while(currentAbsolute<maxAbsoluteValueSquared && iters<maxIter);
  return iters;
}

int main ( int argc, char *argv[] )
{
  int rank,nProc,errorCode;
  double b=2;
  int N=256;
  double minx=-2,maxx=0,miny=-1,maxy=1;
  int width=2048, height=2048;
  double dx=(maxx-minx)/(width-1),dy=(maxy-miny)/(height-1);


  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  //Checking of number of columns is dividable.
  if (width%nProc!=0)
  {
    std::cout << "Number of columns not evenly divisable by the number of processes.";
    MPI_Finalize();
    return 0;
  }

  //Some help variables to simplify code.
  //Picture divided vertically into segments.
  int partitionWidth=width/nProc;
  int partitionSize=height*partitionWidth;
  int resArray[partitionSize];

  //Calculating certain part of the picture for certain processes.
  for(int xIter=rank*partitionWidth;xIter<(rank+1)*partitionWidth; xIter++)
  {
    for(int yIter=0; yIter<height;yIter++)
    {
      //Saving the value calculated by the function calcPix.
      resArray[(xIter%partitionWidth)*height+yIter]=
        calcPix(b,N,((double)(xIter))*dx+minx,((double)(yIter))*dy+miny);
    }
  }

  //Collecting all information
  if (rank==0) {
    FILE *fp;
    fp = fopen("color.txt","w");

    for (int i=0; i<nProc; i++)
    {
      //Fetching calculated arrays from other processes if not "master" process.
      if(i!=0) errorCode=MPI_Recv(&resArray,partitionSize, MPI_INT,i,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int j = 0;j < partitionSize; j++) {
        fprintf(fp, "%hhu", resArray[j]); // printing pixelvalue to file
        fprintf(fp, "\n");
      }
    }
    fclose(fp);
  }
  if(nProc>1&&rank>0)
  {
    MPI_Send(&resArray, partitionSize, MPI_INT,0,1,MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}
