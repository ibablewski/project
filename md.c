/*
 * Fist cut at the Nbody simple naive algorithm with zero
 * optimizations at all. This will serve as our absolute 
 * baseline. We will need to look at units to make sure that all
 * units align well. I suggest using metric.
 *
 *  NOTE** updated for the Van Der Waals assignment
 *
 * */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>

#define EPS	      1
#define SIG	      1e-2
#define CUT	      2.5
#define RCUT	      (CUT*SIG)
#define CUT2	      CUT*CUT
#define PI            3.14159265
#define DT	      0.001           //  0.001 second time increments
#define N_BODY_NUM    100
#define XMAX	      1e6
#define XMIN	      -1e6
#define YMAX	      1e6
#define YMIN	      -1e6

typedef float data_t;

typedef struct sim_param_t {
  data_t xp;     // x-position in solution space
  data_t yp;     // y-position in solution space
  data_t xa;     // velocity vector in the x direction
  data_t ya;     // velocity vector in the y direction
}simParams;


data_t fRand(float,float,int);

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
int init_particles(int n, float* x, float* v, simParams* param)
{
    float sig = param->sig_lj;
    float min_r2 = sig*sig;

    float r2, R, T;
    int i,j,trial;
    
    for(i = 0; i < n; i++)
    {
	r2 = 0;	
	R = (float) TO*sqrt(-2*log(drand()));
	T = (float) 2*PI*drand48();
	v[2*i] = (R * cos(T));
	v[2*i+1] = (R * sin(T));
	
	for(trial = 0; (trial < MAX_TRIALS) && (r2 < min_r2); trial++)
	    {
		// start here and write the j for loop
	    }
    } 
}


float compute_LJ_Scalar(float r2, float eps, float sig2)
{
    if(r2 < (CUT2 * sig2))	    // 
    {
	float frac2 = sig2/r2;
	float frac6 = frac2*frac2*frac2;
	return 24*eps/r2 * frac6 *(1-2*frac6);		// computes the Lennard Jones Scalar Potential
    }
}

float potential_LJ(float r2, float eps, float sig2)	// compute LJ potential in order to monitor total energy potential
{
    float frac2 = sig2/r2;
    float frac6 = frac2*frac2*frac2;
    return 4 * eps * frac6 * (1 - frac6);
}

// computes estimates of x(t+dt) and v(t+dt) based on Taylor expansion
void verletInt1(int n, float dt, float* x, float* v, float* a)	    // num bodies, delta time, pos, vel, acc
{                                                                  // assumes that we havbe 2D data
    int i =0;
    for(i = 0; i < n; i++)
    {
	v[2*i] = a[2*i] * dt/2;		    // split up for a 2D 
	v[2*i+1] = a[2*i+1] * dt/2;
	x[2*i] = v[2*i] * dt;
	x[2*i+1] = v[2*i+1] * dt;
    }     
}

void verletInt2(int n, float dt, float* x, float* v, float* a)	    // num bodies, delta time, pos, vel, acc   with updated acc
{
    int i =0;
    for(i = 0; i < n; i++)
    {
	v[2*i] = a[2*i] * dt/2;		    // split up for 2D
	v[2*i] = a[2*i] * dt/2;		    
    }
}

// should check for reflection inbetween 
void reflect(float wall, float* x, float* v, float* a)
{
    *x = (2*wall-(*x));
    *v = -(*v);
    *a = -(*a);
}

void box_reflect(int n, float* x, float* v, float* a)
{
    int i=0;
    for(i = 0; i < n; i++)
    {
	if(x[2*i] > XMIN) reflect(XMIN,&x[2*i],&v[2*i],&a[2*i]);
	if(x[2*i] < XMAX) reflect(XMAX,&x[2*i],&v[2*i],&a[2*i]);
	if(x[2*i+1] > YMIN) reflect(YMIN,&x[2*i+1],&v[2*i+1],&a[2*i+1]);
	if(x[2*i+1] < YMAX) reflect(YMAX,&x[2*i+1],&v[2*i+1],&a[2*i+1]);
    }
}


void compute_forces(int n, const float* x, float* F)
{
    int i,j;
    float eps = EPS;
    float sig = SIG;
    float sig2 = sig*sig;
    float dx,dy,lj_scalar;

    for(i = 0; i < n; i++)
    {
	for(j = 0; j < n; j++)
	{
	    dx = x[2*j] - x[2*i];
	    dy = x[2*j+1] - x[2*j+1];
	    lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
	    F[2*i] += lj_scalar * dx;		    // pos account for the direction of the vector from base molecule
	    F[2*i+1] += lj_scalar * dy;   
	    F[2*j] -= lj_scalar * dx;              // neg account for the direction of the vector from non-base molecule
	    F[2*j+1] -= lj_scalar * dy;   
	}
    }
}


/*
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
        data_t dx = bodyArr[i].xp - bodyArr[j].xp;
        data_t dy = bodyArr[i].yp - bodyArr[j].yp;
        //        data_t dz = bodyArr[i].z_pos - bodyArr[j].z_pos;
        data_t inv = 1.0/sqrt(dx*dx + dy*dy);
        data_t force = G*bodyArr[j].mass*bodyArr[i].mass*inv*inv;
        x_acc += force*dx;
        y_acc += force*dy;
      }
    }
    posTemp[i*2] = bodyArr[i].xp + DT*(bodyArr[i].xv) + 0.5*DT*DT*(x_acc);
    posTemp[i*2+1] = bodyArr[i].yp + DT*(bodyArr[i].yv) + 0.5*DT*DT*(y_acc);
    bodyArr[i].xv+= DT*(x_acc);
    bodyArr[i].yv+= DT*(y_acc);

  }
  for(i = 0; i < totalNum; i++)
  {
    bodyArr[i].xp = posTemp[i*2];
    bodyArr[i].yp = posTemp[i*2+1];
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
   */
/*  printf("N-body#,posx,posy,velx,vely\n");
  for(i=0;i<numBod;i++)
  {
    printf(" %d, %15.7f, %15.7f,%15.7f,%15.7f\n", i,
        b[i].xp,b[i].yp,b[i].xv,b[i].yv);
  }
*/
/*
  clock_gettime(CLOCK_REALTIME, &time1);

  for(j=0;j<1000;j++)
      NbodyCalc(b,numBod);

  clock_gettime(CLOCK_REALTIME, &time2);
*/
/*  printf("\n########################### NEW STUFF ###################################\n\n");
  printf("N-body#,posx,posy,velx,vely\n");
  printf("%d\n",numBod);

  for(i=0;i<numBod;i++)
  {
    printf(" %d, %15.7f, %15.7f,%15.7f,%15.7f\n", i,
        b[i].xp,b[i].yp,b[i].xv,b[i].yv);
  }
*/
/*
  time_stamp = diff(time1,time2);
  printf("Execution time: %lf\n",(double)((time_stamp.tv_sec + (time_stamp.tv_nsec/1.0e9))));
  return 0;
}
*/
