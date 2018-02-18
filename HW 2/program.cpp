# include <iostream>
# include <mpi.h>
# include <cmath>

using namespace std;

double f(double x)
{
  return 1+x;
}

int main ( int argc, char *argv[] )
{
  int rank,nProc,errorCode;
  int N=99, start=0, end=1;
  double h=1/((double)N +1), uInitial=0, uEnd=0;





  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  int go=1;
  if(rank==0)
  {
      std::cout << "I'm process nr: "<< rank+1<< '\n';
      std::cout << "With points from " <<h*rank*(N+1)/nProc<<" to "<<h*(rank+1)*(N+1)/nProc<< '\n';
      for (int i = rank*(N+1)/nProc;i <(rank+1)*(N+1)/nProc ; i++) {
        std::cout << h*i <<", ";
      }
      std::cout << '\n';

      if(nProc>1)
        MPI_Send(&go, 1, MPI_INT,rank+1,0,MPI_COMM_WORLD);
  }
  else
  {
    MPI_Recv(&go, 1, MPI_INT,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::cout << "I'm process nr: "<< rank+1<< '\n';
    std::cout << "With points from " <<h*rank*(N+1)/nProc<<" to "<<h*(rank+1)*(N+1)/nProc<< '\n';
    for (int i = rank*(N+1)/nProc;i <(rank+1)*(N+1)/nProc ; i++) {
      std::cout << h*i <<", ";
    }
    std::cout <<"\n\n";
    if(rank<nProc-1)
      MPI_Send(&go, 1, MPI_INT,rank+1,0,MPI_COMM_WORLD);

  }





  MPI_Finalize();
  return 0;
}
