/*
 * Fist cut at the Nbody simple naive algorithm with zero
 * optimizations at all. This will serve as our absolute 
 * baseline. We will need to look at units to make sure that all
 * units align well. I suggest using metric.
 *
 * */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>


#define G             6.67384E-11
#define PI            3.14159265
#define DT	      0.001           //  0.001 second time increments
#define N_BODY_NUM    100

typedef float data_t;

typedef struct body {
  data_t mass;    // constant mass
  data_t x_pos;     // x-position in solution space
  data_t y_pos;     // y-position in solution space
  data_t x_vel;     // velocity vector in the x direction
  data_t y_vel;     // velocity vector in the y direction
}Body;

typedef struct forceVec {     // stuct to hold vector data on total Force exerted on single body
  data_t x_vel;             // don't think we need this. I Think we may be able to eliminate
  data_t y_vel;
}Force;

inline data_t invDistance(Body* body1, Body* body2)
{
  return (data_t) sqrt((body1->x_pos-body2->x_pos)*(body1->x_pos-body2->x_pos) +
      (body1->y_pos-body2->y_pos)*(body1->y_pos-body2->y_pos));
}

data_t fRand(float,float,int);

int initBodies(Body* bodies, int Num)
{
  // Initialize all of the bodies
  int i,mult=12827467;
  for(i=0; i< Num; i++)
  {
    bodies[i].mass = fRand(100000.0,-100000.0,i);
    bodies[i].x_pos = fRand(10,-10,i+mult);
    bodies[i].y_pos = fRand(10,-10,i+2*mult);
    bodies[i].x_vel = fRand(10,-10,i+3*mult);
    bodies[i].y_vel = fRand(10,-10,i+4*mult);
    mult+=81722171%192323820820392;
  }

}

data_t fRand(float max, float min, int seed)     // need to iterate over seed key or else same result
{
  srand(seed);
  return (data_t) (min + ((float)(rand()/(float)RAND_MAX)*(float)(max-min)));
}

/*
 *  Want to calculate (initially) the gravitational forces on all N
 *  bodies.
 *  F = G*(m_1*m_2)/r^2  Where G is the gravitational Constant
 *  formally, the force exerted from gravity is operable in the opposite direction of
 *  it's field. Therefore, the above will be nevative in practice.
 *
 */


void NbodyCalc(Body* bodyArr, int totalNum)
{
  int i,j;
  data_t posTemp[totalNum*2];
  for(i=0;i<totalNum;i++)
  {
    data_t x_acc=0;
    data_t y_acc=0;
    //    data_t z_acc=0;
    for(j=0;j<totalNum;j++)
    {
      if(i!=j){
        data_t dx = bodyArr[i].x_pos - bodyArr[j].x_pos;
        data_t dy = bodyArr[i].y_pos - bodyArr[j].y_pos;
        //        data_t dz = bodyArr[i].z_pos - bodyArr[j].z_pos;
        data_t inv = 1.0/sqrt(dx*dx + dy*dy);
        data_t force = G*bodyArr[j].mass*bodyArr[i].mass*inv*inv;
        x_acc += force*dx;
        y_acc += force*dy;
      }
    }
    posTemp[i*2] = bodyArr[i].x_pos + DT*(bodyArr[i].x_vel) + 0.5*DT*DT*(x_acc);
    posTemp[i*2+1] = bodyArr[i].y_pos + DT*(bodyArr[i].y_vel) + 0.5*DT*DT*(y_acc);
    bodyArr[i].x_vel+= DT*(x_acc);
    bodyArr[i].y_vel+= DT*(y_acc);

  }
  for(i = 0; i < totalNum; i++)
  {
    bodyArr[i].x_pos = posTemp[i*2];
    bodyArr[i].y_pos = posTemp[i*2+1];
  }
}

struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

int main()
{

  struct timespec diff(struct timespec start, struct timespec end);
  struct timespec time1, time2;
  struct timespec time_stamp;

  Body b[N_BODY_NUM];
  int i,j;
  int numBod = N_BODY_NUM;
  srand(time(NULL));
  initBodies(b,numBod);

/*  printf("N-body#,posx,posy,velx,vely\n");
  for(i=0;i<numBod;i++)
  {
    printf(" %d, %15.7f, %15.7f,%15.7f,%15.7f\n", i,
        b[i].x_pos,b[i].y_pos,b[i].x_vel,b[i].y_vel);
  }
*/
  clock_gettime(CLOCK_REALTIME, &time1);

  for(j=0;j<1000;j++)
      NbodyCalc(b,numBod);

  clock_gettime(CLOCK_REALTIME, &time2);
/*  printf("\n########################### NEW STUFF ###################################\n\n");
  printf("N-body#,posx,posy,velx,vely\n");
  printf("%d\n",numBod);

  for(i=0;i<numBod;i++)
  {
    printf(" %d, %15.7f, %15.7f,%15.7f,%15.7f\n", i,
        b[i].x_pos,b[i].y_pos,b[i].x_vel,b[i].y_vel);
  }
*/
  time_stamp = diff(time1,time2);
  printf("Execution time: %lf\n",(double)((time_stamp.tv_sec + (time_stamp.tv_nsec/1.0e9))));
  return 0;
}

