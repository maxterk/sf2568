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
  int nParticles=5000;
  double xMin=-1,xMax=1,yMin=-1,yMax=1;
  double minima=sqrt((xMax-xMin)*(yMax-yMin)/nParticles*4/sqrt(3.0));
  double cutOffDistance=minima*3;



  int rank,nProc,errorCode;
  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  if(rank==0)
    std::cout << minima << '\n';

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
      // double theta=randBetween(0,2*3.1415926535);
      // double r=randBetween(0,3);
      particles[i].x=randBetween(xMin,0);
      particles[i].y=randBetween(yMin,0);

      // particles[i].x=r*cos(theta);
      // particles[i].y=r*sin(theta);

      // particles[i].y=r*exp(sin(theta));
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

// if(rank==0)
// {
//   // std::cout << "/* message */" << '\n';
//   // std::cout << indexX<< " and "<< nParticles/xDim << '\n';
// }
  std::nth_element (xVals.begin(), xVals.begin()+(indexX+1)*nParticles/xDim-1, xVals.end());
  std::nth_element (xVals.begin(), xVals.begin()+indexX*nParticles/xDim, xVals.end());
  double xLower=xVals[indexX*nParticles/xDim],
          xUpper=xVals[(indexX+1)*nParticles/xDim-1];

  std::vector<Particle> xSelect;

  for(int i=0; i<nParticles; i++)
  {
    if(xLower <=particles[i].x && particles[i].x<=xUpper)
        yVals.push_back(particles[i].y);

    if(xLower-cutOffDistance<=particles[i].x && particles[i].x<=xUpper+cutOffDistance)
        xSelect.push_back(particles[i]);

  }
  // std::cout << tempSet.size()%yDim << '\n';
  std::nth_element (yVals.begin(), yVals.begin()+(indexY+1)*yVals.size()/yDim-1, yVals.end());
  double yUpper=yVals[(indexY+1)*yVals.size()/yDim-1];

  std::nth_element (yVals.begin(), yVals.begin()+indexY*yVals.size()/yDim, yVals.end());
  double yLower=yVals[indexY*yVals.size()/yDim];


  std::vector<Particle> updateSet;
  std::vector<Particle> complementarySet;

  for(int i=0; i<xSelect.size(); i++)
  {
    // if(rank==0)
    //   std::cout << tempSet[i].y << '\n';
    if(yLower-cutOffDistance<=xSelect[i].y && xSelect[i].y<=yUpper+cutOffDistance)
    {
      if(xLower <=xSelect[i].x && xSelect[i].x<=xUpper)
      {
        if(yLower<=xSelect[i].y && xSelect[i].y<=yUpper)
        {
          updateSet.push_back(xSelect[i]);
        }
        else
        {
          complementarySet.push_back(xSelect[i]);
        }
      }
      else
      {
        complementarySet.push_back(xSelect[i]);
      }
    }
  }

  //Updating particles
  for(int i=0; i < updateSet.size(); i++)
  {
    updateSet[i].x+=randBetween(-0.1,0.1);
    updateSet[i].y+=randBetween(-0.1,0.1);

    updateSet[i].x=min(updateSet[i].x,xMax);
    updateSet[i].x=max(updateSet[i].x,xMin);
    updateSet[i].y=min(updateSet[i].y,yMax);
    updateSet[i].y=max(updateSet[i].y,yMin);
  }


  std::vector<Particle> allParticles;
  int nElems=rank+1;
  int particlesPerProcessor[nProc],
      cumPart[nProc];
  cumPart[0]=0;
  particlesPerProcessor[rank]=nElems;
  Particle particles[nElems];

  for(int i=0; i<nProc; i++)
  {
    MPI_Bcast(&particlesPerProcessor[i], 1, MPI_INT, i, MPI_COMM_WORLD);
    // if(i<nProc-1 && i>0)
    // {
    //   cumPart[i+1]=cumPart[i]+particlesPerProcessor[i];
    // }
    // Particle ps[particlesPerProcessor[i]];
    // if(rank==i)
    // {
    //   for(int j=0; j<particlesPerProcessor[i]; j++)
    //   {
    //     ps[j]=particles[j];
    //   }
    // }
    // MPI_Bcast(&ps,particlesPerProcessor[i]*2, MPI_DOUBLE,i,MPI_COMM_WORLD);
    //
    //
    // allParticles.insert(allParticles.end(),ps,ps+particlesPerProcessor[i]);
  }

if(rank==0)
{
  for (int i =0; i < nProc; i++)
  {
    std::cout << "Particles at processor "<< i+1<<" is "<< particlesPerProcessor[i] << '\n';
  }
}
  //Synchronized printout
//   FILE *fp;
//   int go=1;
//   if(rank!=0)
//   {
//     MPI_Recv(&go, 1, MPI_INT,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     fp = fopen("output.txt","a");
//   }
//   else
//   {
//     fp = fopen("output.txt","w");
//   }
//
//   //Actual printout
//   std::cout << "Greetings from processor: "<<rank+1 << '\n';
//   std::cout << "Currently manages: "<<updateSet.size() << '\n';
//   std::cout << "with complementary "<< complementarySet.size() << '\n';
//   // std::cout << "Xind " <<indexX;
//   // std::cout << " Yind " <<indexY<< '\n';
//   // std::cout << "Xdim " <<xDim;
//   // std::cout << " Ydim " <<yDim<< '\n';
//
//   for(int i=0; i<updateSet.size(); i++)
//   {
//     // std::cout <<"x: "<< updateSet[i].x<<"y: "<< updateSet[i].y<< '\n';
//     fprintf(fp, "%f ", updateSet[i].x);
//     fprintf(fp, "%f ", updateSet[i].y);
//     fprintf(fp, "%f ", (indexX+indexY)%2==1?1.0:0.0);
//   }
//   // if(rank==7)
//   // {
//   //   for(int i=0; i<updateSet.size(); i++)
//   //   {
//   //     // std::cout <<"x: "<< updateSet[i].x<<"y: "<< updateSet[i].y<< '\n';
//   //     fprintf(fp, "%f ", updateSet[i].x);
//   //     fprintf(fp, "%f ", updateSet[i].y);
//   //     fprintf(fp, "%f ", 1.0);
//   //   }
//   //
//   //   for(int i=0; i<complementarySet.size(); i++)
//   //   {
//   //     fprintf(fp, "%f ", complementarySet[i].x);
//   //     fprintf(fp, "%f ", complementarySet[i].y);
//   //     fprintf(fp, "%f ", 0.0);
//   //   }
//   // }
//   // std::cout << "Spacing" << '\n';
//
//   // for(int i=0; i<tempSet.size(); i++)
//   // {
//   //   std::cout <<"x: "<< tempSet[i].x<<"y: "<< tempSet[i].y<< '\n';
//   // }
//   // std::cout << tempSet.size() << '\n';
//   // std::cout << "Y Upper: "<<yUpper << '\n';
//   // std::cout << "Y Lower: "<<yLower << '\n';
//
// //  std::cout << "" << '\n'<<flush;
//   // std::cout << xLower<< " to "<<xUpper << '\n';
//   fflush(stdout);
//
//   //Handing over to next process.
//   fclose(fp);
//   if(nProc>1)
//   {
//     if(rank<nProc-1)
//       MPI_Send(&go, 1, MPI_INT,rank+1,0,MPI_COMM_WORLD);
//   }

  MPI_Finalize();
  return 0;
}
