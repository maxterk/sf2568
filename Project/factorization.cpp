# include <iostream>
# include <cmath>
# include <mpi.h>
# include <algorithm>
# include <vector>

using namespace std;

void factor(int processors, int &a, int &b)
{
  a=(int) sqrt(processors);
  b=processors/a;
  while(a*b!=processors)
  {
    a--;
    b=processors/a;
  }
}

int main ( int argc, char *argv[] )
{
  int rank,nProc,errorCode;
  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  int a,b;
  factor(nProc,a,b);
  int indexY=rank%a;
  int indexX=(rank-indexY)/b;

  std::cout << "Processor: "<<
  rank+1<<"|"<<a<<"|"<<b<<"|"<<
  indexX<<"|"<<indexY << '\n';

  MPI_Finalize();
}
