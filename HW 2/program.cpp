# include <iostream>
# include <mpi.h>
# include <cmath>

using namespace std;

double f(double x)
{
  return sin(5*x);
}
double r(double x)
{
  return -exp(x);
}

int main ( int argc, char *argv[] )
{
  //Problem parameters
  int maxIter=1000000;//1000000
  int N=1000;
  double uFirst=0, uEnd=0, h=1/((double)N+1);

  //cheaper calculation
  double h2=pow(h,2);

  int rank,nProc,errorCode;
  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );


  //Checking how many available processes
  if(nProc>N &&rank==0)
  {
    std::cout << "Error, too many processes" << '\n';
    return 1;
  }

  //Allocating range of indices for the different processes.
  //May be of unequal size.
  int prelWidth=N/nProc;
  int toSpread=N%nProc;
  int start,end=0;
  if(rank<toSpread)
  {
    start=(prelWidth+1)*rank+1;
    end=start+prelWidth;
  }
  else
  {
    start=(prelWidth+1)*(toSpread-1)+prelWidth*(rank-toSpread+1)+2;
    end=start+prelWidth-1;
  }

  //Allocating arrays for iterative solution
  double currentSolution [end-start+3];
  double nextSolution [end-start+3];
  int pointsInArray=end-start+3;

  //Initializing first guess
  for(int i=0; i<pointsInArray; i++)
  {
    currentSolution[i]=1;
  }

  //Making sure initial guess conforms to end values
  if(rank==0) currentSolution[0]=uFirst;
  if(rank==nProc-1) currentSolution[pointsInArray-1]=uEnd;


  for(int iter=0; iter<maxIter; iter++)
  {
      //Calculating next solution
      for(int i=1; i<pointsInArray-1; i++)
      {
          nextSolution[i]=
            (currentSolution[i-1]+currentSolution[i+1]-h2*f((start+i-1)*h))/
              (2.0-h2*r(h*(start+i-1)));
      }

      //As we use old values to calculate new ones
      //we need to keep track of them separately.
      for(int i=1; i<pointsInArray-1; i++)
      {
        currentSolution[i]=nextSolution[i];
      }

      // Communicating values
      if(nProc>1)
      {
        //Even processes' ends overlapping with odds starts
        if(rank%2==0 && rank!=nProc-1) MPI_Send(&currentSolution[pointsInArray-2], 1, MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
        if(rank%2==1) MPI_Recv(&currentSolution[0], 1, MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(rank%2==0 && rank!=nProc-1) MPI_Recv(&currentSolution[pointsInArray-1], 1, MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(rank%2==1) MPI_Send(&currentSolution[1], 1, MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);

        //odd processes' ends overlapping with evens starts
        //Left part of overlap
        if(rank%2==1 && rank!=nProc-1) MPI_Send(&currentSolution[pointsInArray-2], 1, MPI_DOUBLE,rank+1,2,MPI_COMM_WORLD);
        if(rank%2==0 && rank!=0) MPI_Recv(&currentSolution[0], 1, MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // //Right part of overlap
        if(rank%2==1 && rank!=nProc-1) MPI_Recv(&currentSolution[pointsInArray-1], 1, MPI_DOUBLE,rank+1,3,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(rank%2==0 && rank!=0) MPI_Send(&currentSolution[1], 1, MPI_DOUBLE,rank-1,3,MPI_COMM_WORLD);
      }
  }

  //Dummy variable for synchronization purposes.
  int go=1;

  //Writing results to file
  FILE *fp;
  if(rank==0)
  {
    fp = fopen("output.txt","w");
  }
  else
  {
    fp = fopen("output.txt","a");
    MPI_Recv(&go, 1, MPI_INT,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  for(int i=rank==0?0:1; i<pointsInArray-1; i++)
  {
    fprintf(fp, "%f ", currentSolution[i]);
    fprintf(fp, "\n");
  }
  if(rank==nProc-1)
    fprintf(fp, "%f ", currentSolution[pointsInArray-1]);
  fclose(fp);
  if(rank==0)
  {
      if(nProc>1)
        MPI_Send(&go, 1, MPI_INT,rank+1,0,MPI_COMM_WORLD);
  }
  else
  {
    if(rank<nProc-1)
      MPI_Send(&go, 1, MPI_INT,rank+1,0,MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
