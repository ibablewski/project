/*
 * Fist cut at the Nbody simple naive algorithm with zero
 * optimizations at all. This will serve as our absolute 
 * baseline.
 *
 * */


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>

// First we need to write up a skeleton algorithm

#define G             6.67384E-11
#define TIME_STEP     0.500           //  0.5 second time increments

typedef struct body {
  const float mass;    // constant mass
  float x_pos;     // x-position in solution space
  float y_pos;     // y-position in solution space
  float u_vec;     // u-vecor = velocity vector in the x direction
  float v_vec;     // v-vector = velocity vector in the y direction

}Body;

typedef struct forceVec {     // stuct to hold vector data on total Force exerted on single body
  float u_vec;
  float v_vec;
}Force;

int initBodies(float* bodies, int Num)
{
  // Initialize all of the bodies
  
}

float fRand(float max, float min, int seed)     // need to iterate over seed key or else same result
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

// -> is shorthand for (*foo).x=46;

int updateBody(Body* pointForce,Force* totalForce)
{
  // recalculate the velocity vectors
  pointForce->u_vec += totalForce->u_vec;
  pointForce->v_vec += totalForce->v_vec;

  // update to new position from new vectors
  pointForce->x_pos +=(TIME_STEP*(pointForce->u_vec));
  pointForce->y_pos +=(TIME_STEP*(pointForce->v_vec));

  return 0;
}


int main()
{
  Body bod;
  float x,y,u,v,r,q;
  u = 10;
  v = 9;
  x = 1;
  y = 2;
  bod.x_pos = x;
  bod.y_pos = y;
  bod.u_vec = u;
  bod.v_vec = v;

  printf("%f %f %f %f \n", bod.u_vec, bod.v_vec, bod.x_pos ,bod.y_pos);

  Force fc;
  fc.u_vec = 7;
  fc.v_vec = 5;

  int err = updateBody(&bod,&fc);
  if(err)
     printf("%d\n", err);
  
  printf("%f %f %f %f \n", bod.u_vec, bod.v_vec, bod.x_pos ,bod.y_pos);
  return 0;
}

