# include <iostream>
# include <mpi.h>
# include <cmath>
# include <algorithm>
# include <vector>

using namespace std;

struct Particle{
  double a,b;
};

int randIntBetween(int a, int b)
{
  return rand()%(b-a)+a;
}

double randBetween(double fMin, double fMax)
{
    double f = (double) rand() / (1.0 + RAND_MAX);
    return fMin + f * (fMax - fMin);
}

int main ( int argc, char *argv[] )
{
    int rank,nProc,errorCode;
  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  // if(rank==0)
  // {
  //   Particle p[2];
  //   p[0].a=1;
  //   p[0].b=2;
  //   p[1].a=3;
  //   p[1].b=4;
  //   MPI_Send(&p,4,MPI_DOUBLE,1,1, MPI_COMM_WORLD);
  // }
  // if(rank==1)
  // {
  //   Particle p[2];
  //   MPI_Recv(&p,4,MPI_DOUBLE, 0,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   std::cout << p[0].a << '\n'<<p[0].b<<'\n'<<p[1].a<<'\n'<<p[1].b<<'\n';
  // }

  std::vector<Particle> allParticles;
  int nElems=rank+1;
  //std::cout << nElems << '\n';
  int particlesPerProcessor[nProc],
      cumPart[nProc];
  cumPart[0]=0;
  particlesPerProcessor[rank]=nElems;
  Particle particles[nElems];
  for(int i=0; i< nElems; i++)
  {
    particles[i].a=rank+1;
    particles[i].b=rank+1;
  }


  for(int i=0; i<nProc; i++)
  {
    MPI_Bcast(&particlesPerProcessor[i], 1, MPI_INT, i, MPI_COMM_WORLD);
    if(i<nProc-1 && i>0)
    {
      cumPart[i+1]=cumPart[i]+particlesPerProcessor[i];
    }
    Particle ps[particlesPerProcessor[i]];
    if(rank==i)
    {
      for(int j=0; j<particlesPerProcessor[i]; j++)
      {
        ps[j]=particles[j];
      }
    }
    MPI_Bcast(&ps,particlesPerProcessor[i]*2, MPI_DOUBLE,i,MPI_COMM_WORLD);


    allParticles.insert(allParticles.end(),ps,ps+particlesPerProcessor[i]);
  }

  int go=1;
  if(rank==0)
  {
    //std::cout << "This is known by process: "<<rank+1 << '\n';
    // std::cout << particlesPerProcessor[rank] << '\n';
    for (int i = 0; i < allParticles.size(); i++) {
      std::cout << allParticles[i].a;
      // if(i<allParticles.size()-1)
      // {
      //   if(allParticles[i].a<allParticles[i+1].a)
      //   {
      //     std::cout << '\n';
      //   }
    }
    std::cout << '\n';
  }
  else
  {
    MPI_Recv(&go, 1, MPI_INT,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // std::cout << "This is known by process: "<<rank+1 << '\n';
    // std::cout << particlesPerProcessor[rank] << '\n';
    for (int i = 0; i < allParticles.size(); i++) {
      std::cout << allParticles[i].a;
      // if(i<allParticles.size()-1)
      // {
      //   if(allParticles[i].a<allParticles[i+1].a)
      //   {
      //     std::cout << '\n';
      //   }

    }
    std::cout << '\n';
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
