# include <iostream>
# include <mpi.h>
# include <cmath>
# include <algorithm>
# include <vector>

using namespace std;

struct Particle{
  double x,y;
};

double randBetween(double fMin, double fMax)
{
    double f = (double) rand() / (1.0 + RAND_MAX);
    return fMin + f * (fMax - fMin);
}
int randIntBetween(double a, double b)
{
  return fmod(rand(),b-a)+a;
}
double median(double a, double b, double c)
{
  double maxVal=max(max(a,b),c);
  double minVal=min(min(a,b),c);
  if(a!=maxVal && a!=minVal) return a;
  if(b!=maxVal && b!=minVal) return b;
  return c;
}

void factor(int nProc, int &a, int &b)
{
  a=(int) sqrt(nProc);
  b=nProc/a;
  while(a*b!=nProc)
  {
    a--;
    b=nProc/a;
  }
}

int main ( int argc, char *argv[] )
{
  //Problem parameters
  int nParticles=1000;
  int cutOffDistance=0.01;
  double xMin=-1,xMax=1,yMin=-1,yMax=1;

  int rank,nProc,errorCode;
  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  //Adjusting amount of particles
  nParticles=nParticles-nParticles%nProc;

  //Setting up grid
  int yDim,xDim;
  factor(nProc,yDim,xDim);
  int indexY=rank%yDim;
  int indexX=rank/yDim;
  if(rank==0)
  {
    std::cout << "      Grid size" << '\n';
    std::cout << "######################" << '\n';
    std::cout << "X-dimension: "<< xDim << '\n';
    std::cout << "Y-dimension: "<< yDim<< '\n';
    std::cout << "######################" << '\n';
  }

  //Setting seed
  //srand(rank);

  //Array for particles
  Particle particles[nParticles];

  //Creating a proper amount of particles on first processor.
  if(rank==0)
  {
    for (int i = 0; i < nParticles; i++)
    {
      particles[i].x=randBetween(xMin,xMax);
      particles[i].y=randBetween(yMin,yMax);
    }
  }
  MPI_Bcast(&particles, nParticles*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::vector<double> xVals;
  std::vector<double> yVals;


  for(int i=0; i<nParticles; i++)
    xVals.push_back(particles[i].x);


  // if (rank==0)
  // {
  //   //std::nth_element (xVals.begin(), xVals.begin(), xVals.end());
  //   for(int i=0; i<nParticles; i++)
  //   {
  //     std::nth_element (xVals.begin(), xVals.begin()+i, xVals.end());
  //     std::cout << xVals[i] << '\n';
  //   }
  // }
  //std::nth_element (xVals.begin(), xVals.begin()+, xVals.end());
  // std::nth_element (xVals.begin(), xVals.begin()+i, xVals.end());

if(rank==0)
{
  // std::cout << "/* message */" << '\n';
  // std::cout << indexX<< " and "<< nParticles/xDim << '\n';
}
  std::nth_element (xVals.begin(), xVals.begin()+(indexX+1)*nParticles/xDim-1, xVals.end());
  std::nth_element (xVals.begin(), xVals.begin()+indexX*nParticles/xDim, xVals.end());
  double xLower=xVals[indexX*nParticles/xDim],
          xUpper=xVals[(indexX+1)*nParticles/xDim-1];

  std::vector<Particle> tempSet;
  for(int i=0; i<nParticles; i++)
  {
    if(xLower <=particles[i].x && particles[i].x<=xUpper)
    {
      tempSet.push_back(particles[i]);
      yVals.push_back(particles[i].y);
    }
  }
  // std::cout << tempSet.size()%yDim << '\n';
  std::nth_element (yVals.begin(), yVals.begin()+(indexY+1)*tempSet.size()/yDim-1, yVals.end());
  double yUpper=yVals[(indexY+1)*tempSet.size()/yDim-1];

  std::nth_element (yVals.begin(), yVals.begin()+indexY*tempSet.size()/yDim, yVals.end());
  double yLower=yVals[indexY*tempSet.size()/yDim];


  std::vector<Particle> updateSet;
  // if(rank==0)
  //   std::cout << "Y-values" << '\n';
  for(int i=0; i<tempSet.size(); i++)
  {
    // if(rank==0)
    //   std::cout << tempSet[i].y << '\n';
    if(yLower <=tempSet[i].y && tempSet[i].y<=yUpper)
    {
      updateSet.push_back(tempSet[i]);
    }
  }
  std::vector<Particle> complementarySet;


  //Synchronized printout
  FILE *fp;
  int go=1;
  if(rank!=0)
  {
    MPI_Recv(&go, 1, MPI_INT,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    fp = fopen("output.txt","a");
  }
  else
  {
    fp = fopen("output.txt","w");
  }

  //Actual printout
  std::cout << "Greetings from processor: "<<rank+1 << '\n';
  std::cout << "Currently manages: "<<updateSet.size() << '\n';
  // std::cout << "Xind " <<indexX;
  // std::cout << " Yind " <<indexY<< '\n';
  // std::cout << "Xdim " <<xDim;
  // std::cout << " Ydim " <<yDim<< '\n';
  for(int i=0; i<updateSet.size(); i++)
  {
    // std::cout <<"x: "<< updateSet[i].x<<"y: "<< updateSet[i].y<< '\n';
    fprintf(fp, "%f ", updateSet[i].x);
    fprintf(fp, "%f ", updateSet[i].y);
    fprintf(fp, "%f ", (indexX+indexY)%2==1?1.0:0.0);
  }
  // std::cout << "Spacing" << '\n';

  // for(int i=0; i<tempSet.size(); i++)
  // {
  //   std::cout <<"x: "<< tempSet[i].x<<"y: "<< tempSet[i].y<< '\n';
  // }
  // std::cout << tempSet.size() << '\n';
  // std::cout << "Y Upper: "<<yUpper << '\n';
  // std::cout << "Y Lower: "<<yLower << '\n';

//  std::cout << "" << '\n'<<flush;
  // std::cout << xLower<< " to "<<xUpper << '\n';
  fflush(stdout);

  //Handing over to next process.
  fclose(fp);
  if(nProc>1)
  {
    if(rank<nProc-1)
      MPI_Send(&go, 1, MPI_INT,rank+1,0,MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
