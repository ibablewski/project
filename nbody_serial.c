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
#define TIME_STEP     0.500           //  0.5 second time increments

typedef float data_t;

typedef struct body {
  const data_t mass;    // constant mass
  data_t x_pos;     // x-position in solution space
  data_t y_pos;     // y-position in solution space
  data_t u_vec;     // u-vecor = velocity vector in the x direction
  data_t v_vec;     // v-vector = velocity vector in the y direction

}Body;

typedef struct forceVec {     // stuct to hold vector data on total Force exerted on single body
  data_t u_vec;
  data_t v_vec;
}Force;

inline data_t distance(Body* body1, Body* body2)
{
  return (data_t) sqrt((body1->x_pos-body2->x_pos)*(body1->x_pos-body2->x_pos) + 
      (body1->y_pos-body2->y_pos)*(body1->y_pos-body2->y_pos));
}
int initBodies(Body* bodies, int Num)
{
  // Initialize all of the bodies

}

data_t fRand(data_t max, data_t min, int seed)     // need to iterate over seed key or else same result
{
  srand(seed);
  return min + ((rand()/RAND_MAX)*(max-min));
}

/*
 *  Want to calculate (initially) the gravitational forces on all N
 *  bodies.
 *  F = G*(m_1*m_2)/r^2  Where G is the gravitational Constant
 *  formally, the force exerted from gravity is operable in the opposite direction of
 *  it's field. Therefore, the above will be nevative in practice.
 * 
 * */


// where body one is point mass, body two is body acting on it
Force calcForce(Body* body1, Body* body2)
{
  Force tempForce;
  data_t dist = distance(body1,body2);
  data_t angle = (data_t) asin(abs(body1->y_pos-body2->y_pos)/dist);
  data_t grav_force = ((-1)* G * (body1->mass * body2->mass))/dist;
  tempForce.u_vec += grav_force * cos(angle);
  tempForce.v_vec += grav_force * sin(angle);

  return tempForce;
}

// calcs total force on one body from all others
// iterating from 0 to N-1
Force calcTotalForce(Body* bodyArr,int bodyNum, int totalNum)
{
  int i,j;
  Force totalForce;
  Force tempForce;
  for(i=0;i<totalNum;i++)
  {
    if(i!=bodyNum){
      tempForce = calcForce(&bodyArr[bodyNum],&bodyArr[i]);     // need to declare this
      totalForce.u_vec += tempForce.u_vec;
      totalForce.v_vec += tempForce.v_vec;
    }
  }
  return totalForce;
}

// don't need this function any more

void updateBody(Body* pointForce)
{
  // recalculate the velocity vectors
  //pointForce->u_vec += totalForce->u_vec;
  //pointForce->v_vec += totalForce->v_vec;

  // update to new position from new vectors
  pointForce->x_pos +=(TIME_STEP*(pointForce->u_vec));
  pointForce->y_pos +=(TIME_STEP*(pointForce->v_vec));

}

void NbodyCalc(Body* bodyArr, int totalNum)
{
  Force tempForce;
  int i = 0;
  for(i=0; i < totalNum;i++)
  {
    tempForce = calcTotalForce(bodyArr, i , totalNum);
    bodyArr[i].u_vec += tempForce.u_vec;
    bodyArr[i].v_vec += tempForce.v_vec;
  }

  for(i=0; i< totalNum; i++)
  {
    bodyArr[i].x_pos +=(TIME_STEP*(bodyArr[i].u_vec));
    bodyArr[i].y_pos +=(TIME_STEP*(bodyArr[i].v_vec));
  }
}

int main()
{
  Body bod,bd;
  data_t x,y,u,v,r,q;
  u = 10;
  v = 9;
  x = 1;
  y = 2;
  bod.x_pos = x;
  bod.y_pos = y;
  bod.u_vec = u;
  bod.v_vec = v;

  bd.x_pos = 8;
  bd.y_pos = 3;
  bd.u_vec = 3;
  bd.v_vec = 2;

  printf("%f %f %f %f \n", bod.u_vec, bod.v_vec, bod.x_pos ,bod.y_pos);

  Force fc;
  fc.u_vec = 7;
  fc.v_vec = 5;

  data_t dis = distance(&bod, &bd);

  updateBody(&bod);

  printf("%f %f %f %f %f\n", bod.u_vec, bod.v_vec, bod.x_pos ,bod.y_pos, dis);
  return 0;
}

