# include <iostream>
# include <mpi.h>
# include <cmath>
# include <algorithm>
# include <vector>

using namespace std;

struct Particle{
  double x,y,
          vx,vy;
};

double randBetween(double fMin, double fMax)
{
    double f = (double) rand() / (1.0 + RAND_MAX);
    return fMin + f * (fMax - fMin);
}
int randIntBetween(int a, int b)
{
  return rand()%(b-a)+a;
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
  int nParticles=10000;
  int cutOffDistance=0.01;
  double xMin=-1,xMax=1,yMin=-1,yMax=1;

  int rank,nProc,errorCode;
  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );


  //Setting up grid
  int a,b;
  factor(nProc,a,b);
  int indexY=rank%a;
  int indexX=(rank-indexY)/b;
  if(rank==0)
  {
    std::cout << "      Grid size" << '\n';
    std::cout << "######################" << '\n';
    std::cout << "X-dimension: "<< a << '\n';
    std::cout << "Y-dimension: "<< b<< '\n';
    std::cout << "######################" << '\n';
  }

  //Setting seed
  srand(rank);

  //Creating a proper amount of particles on a processor.
  Particle particles[nParticles];
  for (int i = 0; i < nParticles; i++)
  {
    particles[i].x=randBetween(xMin,xMax);
    particles[i].y=randBetween(yMin,yMax);
  }

  double xMedian=median(
    particles[randIntBetween(0,nParticles)].x,
    particles[randIntBetween(0,nParticles)].x,
    particles[randIntBetween(0,nParticles)].x);

  double yMedian=median(
    particles[randIntBetween(0,nParticles)].y,
    particles[randIntBetween(0,nParticles)].y,
    particles[randIntBetween(0,nParticles)].y);

    double counter=0;
    std::vector<Particle> lxly;
    std::vector<Particle> lxsy;
    std::vector<Particle> sxly;
    std::vector<Particle> sxsy;

    for (int i = 0; i < nParticles; i++)
    {
      if(particles[i].x>xMedian)
      {
        if(particles[i].y>yMedian)
        {
          lxly.push_back(particles[i]);
        }
        else
        {
          lxsy.push_back(particles[i]);
        }
      }
      else
      {
        if(particles[i].y>yMedian)
        {
          sxly.push_back(particles[i]);
        }
        else
        {
          sxsy.push_back(particles[i]);
        }
      }
    }


  int go=1;

  //Writing results to file

  if(rank==0)
  {
    std::cout << lxly.size() << ", "
      << lxsy.size() <<  ", "
      << sxly.size() <<  ", "
      << sxsy.size()<< '\n';
  }
  else
  {
    std::cout << lxly.size() << ", "
      << lxsy.size() <<  ", "
      << sxly.size() <<  ", "
      << sxsy.size()<< '\n';
    MPI_Recv(&go, 1, MPI_INT,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
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
