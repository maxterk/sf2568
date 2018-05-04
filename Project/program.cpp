# include <iostream>
# include <mpi.h>
# include <cmath>
# include <algorithm>
# include <vector>
# include <string>

using namespace std;

struct Particle{
  double x,y,vx,vy;
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
double potentialDerivative(double dist, double minima,double decay)
{
  double eps=4;
  if(dist==minima)
    return 0;

  return 12*eps*(pow(minima/dist,6)/pow(dist,2)-pow(minima/dist,12)/pow(dist,2));
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
  double xMin=-1,xMax=1,yMin=-1,yMax=1;
  double minima=sqrt((xMax-xMin)*(yMax-yMin)/nParticles*4/sqrt(3.0));
    minima=0.1;
  double cutOffDistance=minima*5;
  double deltaT=0.01, tMax=100;

  int rank,nProc,errorCode;
  //Initializing MPI and getting world size and rank of process.
  errorCode = MPI_Init ( &argc, &argv );
  errorCode = MPI_Comm_size ( MPI_COMM_WORLD, &nProc );
  errorCode = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  // Setting up some variables
  int particlesPerProcessor[nProc], cumPart[nProc];
  cumPart[0]=0;
  std::vector<Particle> updateSet;
  std::vector<Particle> complementarySet;
  std::vector<double> xVals;
  std::vector<double> yVals;
  std::vector<Particle> xSelect;
  double xLower,xUpper, yLower, yUpper;

  //Array for particles
  Particle particles[nParticles];

  // if(rank==0)
  //   std::cout<<"Calculated minima: " << minima << '\n';

  //Adjusting amount of particlesparticles
  nParticles=nParticles-nParticles%nProc;

  //Setting up grid
  int yDim,xDim;
  factor(nProc,yDim,xDim);
  int indexY=rank%yDim;
  int indexX=rank/yDim;
  // if(rank==0)
  // {
  //   std::cout << "      Grid size" << '\n';
  //   std::cout << "######################" << '\n';
  //   std::cout << "X-dimension: "<< xDim << '\n';
  //   std::cout << "Y-dimension: "<< yDim<< '\n';
  //   std::cout << "######################" << '\n';
  //   std::cout << "N particles: "<<nParticles << '\n';
  // }

  //Setting seed
  //srand(rank);

  //For timekeeping
  std::vector<double> iterationTimes;

  //Creating a proper amount of particles on first processor.
  if(rank==0)
  {
    for (int i = 0; i < nParticles; i++)
    {
      particles[i].x=randBetween(xMin,xMax);
      particles[i].y=randBetween(yMin,yMax);

      particles[i].vx=0;
      particles[i].vy=0;
    }
  }
  MPI_Bcast(&particles, nParticles*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);


//Main update loop
for(double t=0; t<tMax+1; t++)
{
  MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
  double start = MPI_Wtime(),time,globalTime;
  // MPI_Reduce(&start, &globalStart, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  //Clearing all sets related to data structure
  xVals.clear();
  yVals.clear();
  updateSet.clear();
  complementarySet.clear();
  xSelect.clear();

  //Adding x-values for splitting
  for(int i=0; i<nParticles; i++)
    xVals.push_back(particles[i].x);

  //Finding split values
  std::nth_element (xVals.begin(), xVals.begin()+(indexX+1)*nParticles/xDim-1, xVals.end());
  std::nth_element (xVals.begin(), xVals.begin()+indexX*nParticles/xDim, xVals.end());
  xLower=xVals[indexX*nParticles/xDim];
  xUpper=xVals[(indexX+1)*nParticles/xDim-1];


  //Adding y-values for y-splitting
  for(int i=0; i<nParticles; i++)
  {
    if(xLower <=particles[i].x && particles[i].x<=xUpper)
        yVals.push_back(particles[i].y);

    if(xLower-cutOffDistance<=particles[i].x && particles[i].x<=xUpper+cutOffDistance)
        xSelect.push_back(particles[i]);

  }

  //Splitting on y-values
  std::nth_element (yVals.begin(), yVals.begin()+(indexY+1)*yVals.size()/yDim-1, yVals.end());
  yUpper=yVals[(indexY+1)*yVals.size()/yDim-1];
  std::nth_element (yVals.begin(), yVals.begin()+indexY*yVals.size()/yDim, yVals.end());
  yLower=yVals[indexY*yVals.size()/yDim];


  //Dividing into update and complementary set
  for(int i=0; i<xSelect.size(); i++)
  {
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
    for(int j=i+1; j<updateSet.size(); j++)
    {
        double dist=sqrt(pow(updateSet[i].x-updateSet[j].x,2)+
          pow(updateSet[i].y-updateSet[j].y,2));

        if(dist<cutOffDistance)
        {
          double dudr=potentialDerivative(dist, minima, t/tMax);
          double
            xixj=updateSet[i].x-updateSet[j].x,
            yiyj=updateSet[i].y-updateSet[j].y,
            fx=xixj*dudr,
            fy=yiyj*dudr;

            updateSet[i].vx+=fx*deltaT;
            updateSet[i].vy+=fy*deltaT;
            updateSet[j].vx-=fx*deltaT;
            updateSet[j].vy-=fy*deltaT;

        }


    }
    for(int j=0; j<complementarySet.size(); j++)
    {
      double dist=sqrt(pow(updateSet[i].x-complementarySet[j].x,2)+
        pow(updateSet[i].y-complementarySet[j].y,2));

      if(dist<cutOffDistance)
      {
        double dudr=potentialDerivative(dist, minima, t/tMax);
        double
          xixj=updateSet[i].x-complementarySet[j].x,
          yiyj=updateSet[i].y-complementarySet[j].y,
          fx=xixj*dudr,
          fy=yiyj*dudr;
          updateSet[i].vx+=fx*deltaT;
          updateSet[i].vy+=fy*deltaT;

      }
    }
  }

  //Updating positions based on speed.
  for(int i=0; i < updateSet.size(); i++)
  {
    updateSet[i].x+=updateSet[i].vx*deltaT;
    updateSet[i].y+=updateSet[i].vy*deltaT;
  }


  //######################################################################
  //Measuring runtime of one iteration in parallel.
  // time=MPI_Wtime()-start;
  // MPI_Reduce(&time, &globalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  // if(rank==0)
  //   iterationTimes.push_back(globalTime);
  //######################################################################

  //synchronizing particles over all processors.
  particlesPerProcessor[rank]=updateSet.size();
  for(int i=0; i<nProc; i++)
  {
    MPI_Bcast(&particlesPerProcessor[i], 1, MPI_INT, i, MPI_COMM_WORLD);
    if(i<nProc-1)
    {
      cumPart[i+1]=cumPart[i]+particlesPerProcessor[i];
    }
    Particle ps[particlesPerProcessor[i]];
    if(rank==i)
    {
      for(int j=0; j<particlesPerProcessor[i]; j++)
      {
        ps[j]=updateSet[j];
      }
    }
    MPI_Bcast(&ps,particlesPerProcessor[i]*4, MPI_DOUBLE,i,MPI_COMM_WORLD);

    for(int j=0; j<particlesPerProcessor[i]; j++)
    {
      particles[cumPart[i]+j]=ps[j];
    }
  }


  //######################################################################
  //Measuring runtime of one iteration in parallel.
  time=MPI_Wtime()-start;
  MPI_Reduce(&time, &globalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank==0)
    iterationTimes.push_back(globalTime);
  //######################################################################
}



  // Synchronized printout
  // FILE *fp;
  // int go=1;
  // if(rank!=0)
  // {
  //   MPI_Recv(&go, 1, MPI_INT,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   // fp = fopen("output.txt","a");
  //   fp = fopen("times.txt","a");
  // }
  // else
  // {
  //   // fp = fopen("output.txt","w");
  //   fp = fopen("times.txt","a");
  // }

  //Actual printout
  // std::cout << "Greetings from processor: "<<rank+1 << '\n';
  // std::cout << "Currently manages: "<<updateSet.size() << '\n';
  // std::cout << "with complementary "<< complementarySet.size() << '\n';
  // std::cout << "Xind " <<indexX;
  // std::cout << " Yind " <<indexY<< '\n';
  // std::cout << "Xdim " <<xDim;
  // std::cout << " Ydim " <<yDim<< '\n';

if(rank==0)
{
  double avg=0;
  for(int i=5; i< iterationTimes.size(); i++)
  {
    avg+=iterationTimes[i];
  }
  std::cout << "Average time: "<<avg/(iterationTimes.size()-5) << '\n';

}

  // for(int i=0; i<updateSet.size(); i++)
  // {
  //   // std::cout <<"x: "<< updateSet[i].x<<"y: "<< updateSet[i].y<< '\n';
  //   fprintf(fp, "%f ", updateSet[i].x);
  //   fprintf(fp, "%f ", updateSet[i].y);
  //   fprintf(fp, "%f ", (indexX+indexY)%2==1?1.0:0.0);
  // }
  // if(rank==7)
  // {
  //   for(int i=0; i<updateSet.size(); i++)
  //   {
  //     // std::cout <<"x: "<< updateSet[i].x<<"y: "<< updateSet[i].y<< '\n';
  //     fprintf(fp, "%f ", updateSet[i].x);
  //     fprintf(fp, "%f ", updateSet[i].y);
  //     fprintf(fp, "%f ", 1.0);
  //   }
  //
  //   for(int i=0; i<complementarySet.size(); i++)
  //   {
  //     fprintf(fp, "%f ", complementarySet[i].x);
  //     fprintf(fp, "%f ", complementarySet[i].y);
  //     fprintf(fp, "%f ", 0.0);
  //   }
  // }200
  // std::cout << "Spacing" << '\n';particles

  // for(int i=0; i<tempSet.size(); i++)
  // {
  //   std::cout <<"x: "<< tempSet[i].x<<"y: "<< tempSet[i].y<< '\n';
  // }
  // std::cout << tempSet.size() << '\n';
  // std::cout << "Y Upper: "<<yUpper << '\n';
  // std::cout << "Y Lower: "<<yLower << '\n';

//  std::cout << "" << '\n'<<flush;
  // std::cout << xLower<< " to "<<xUpper << '\n';
  // fflush(stdout);
  //
  // //Handing over to next process.
  // fclose(fp);
  // if(nProc>1)
  // {
  //   if(rank<nProc-1)
  //     MPI_Send(&go, 1, MPI_INT,rank+1,0,MPI_COMM_WORLD);
  // }

  MPI_Finalize();
  return 0;
}
