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

// First we need to write up a skeleton algorithm

#define G     6.67384E-11

struct body {
  const float* mass;    // constant mass
  float* x_pos;     // x-position in solution space
  float* y_pos;     // y-position in solution space
  float* u_vec;     // u-vecor = velocity vector in the x direction
  float* v_vec;     // v-vector = velocity vector in the y direction

};

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


