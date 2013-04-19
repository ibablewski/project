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
#include <string.h>

#define EPS	      1
#define SIG	      1e-2
#define CUT	      2.5
#define RCUT	      (CUT*SIG)
#define CUT2	      CUT*CUT
#define PI            3.14159265
#define DT	      0.001           //  0.001 second time increments
#define N_BODY_NUM    100
#define XMAX	      0.5
#define XMIN	      -0.5
#define YMAX	      0.5
#define YMIN	      -0.5
#define T0          1
#define MAX_TRIALS  100
#define ITERS       100
#define BOX_SIZE    1.0

typedef float data_t;

typedef struct sim_param_t {
  int npart;
  float dt;
  float eps_lj;
  float sig_lj;
} params;

typedef struct molecule {
  float* x;
  float* v;
  float* a;
  float* F;
}mols;

data_t fRand(float,float,int);

data_t fRand(float max, float min, int seed)     // need to iterate over seed key or else same result
{
  srand(seed);
  return (data_t) (min + ((float)(rand()/(float)RAND_MAX)*(float)(max-min)));
}

/*
 *
 *
 *
 *
 *
 *
 */

int init_particles(int n, float* x, float* v, params* param)
{
  float sig = param->sig_lj;
  float min_r2 = sig*sig;

  float r2,dx,dy;
  int i,j,trial;

  for(i = 0; i < n; i++)
  {
    r2 = 0;	

    /* Choose new point via rejection sampling */
    for(trial = 0; (trial < MAX_TRIALS) && (r2 < min_r2); trial++)
    {
      x[2*i] = (float) (BOX_SIZE*drand48()) - BOX_SIZE/2.0;
      x[2*i+1] = (float) (BOX_SIZE*drand48()) - BOX_SIZE/2.0;

      for(j=0; j < i; j++)
      {
        dx = x[2*i] - x[2*j];
        dy = x[2*i+1] - x[2*j+1];
        r2 = dx*dx + dy*dy;
        //printf("Sample%d:%d %f %f %f %f\n",i,j,min_r2,r2,dx,dy);
        if(r2 < min_r2)
          break;
      }
    }
    /* If it takes too many trials, bail and declare number set up */
    if(i > 0 && r2 < min_r2)
      return i;
  }
  return n;
}

void init_particles_va(int n, float* v,float* a, params* param)
{
  float R,T;
  int i;
  
  for(i=0; i < n; i++)
  {
    R = T0 * sqrt(-2.0 * log(drand48()));
    T = 2 * PI * drand48();
    v[2*i] = (R * cos(T));
    v[2*i+1] = (R * sin(T));
    //    printf("SampleVel%d %f %f\n",i,v[2*i],v[2*i+1]);
    a[2*i] = (R * cos(T))/param->dt;
    a[2*i+1] = (R * sin(T))/param->dt;
  }
}


float compute_LJ_Scalar(float r2, float eps, float sig2)
{
  if(r2 < (CUT2 * sig2))	    // 
  {
    float frac2 = sig2/r2;
    float frac6 = frac2*frac2*frac2;
    return 24.0*eps/r2 * frac6 *(1.0-2.0*frac6);		// computes the Lennard Jones Scalar Potential
  }
}

float potential_LJ(float r2, float eps, float sig2)	// compute LJ potential in order to monitor total energy potential
{
  float frac2 = sig2/r2;
  float frac6 = frac2*frac2*frac2;
  return 4.0 * eps * frac6 * (1.0 - frac6);
}

// computes estimates of x(t+dt) and v(t+dt) based on Taylor expansion
void verletInt1(int n, float dt, float* x, float* v, float* a)	    // num bodies, delta time, pos, vel, acc
{                                                                  // assumes that we havbe 2D data
  int i =0;
  for(i = 0; i < n; i++)
  {
    v[2*i] = a[2*i] * dt/2.0;		    // split up for a 2D 
    v[2*i+1] = a[2*i+1] * dt/2.0;
    x[2*i] = v[2*i] * dt;
    x[2*i+1] = v[2*i+1] * dt;
  }
}

void verletInt2(int n, float dt, float* x, float* v, float* a)	    // num bodies, delta time, pos, vel, acc   with updated acc
{
  int v0,v1,i =0;
  for(i = 0; i < n; i++)
  {
    v0 = v[2*i];
    v1 = v[2*i+1];
    v[2*i] = a[2*i] * dt/2.0;		    // split up for 2D
    v[2*i] = a[2*i] * dt/2.0;
    a[2*i] += (v[2*i]-v0)/dt;
    a[2*i+1] += (v[2*i+1]-v1)/dt;
  }
}

void computeAcc(int n, float* F, float* a)        // mass assumed to my constant and 1
{
  int i;
  for(i=0; i< n; i++)
  {
    a[2*i] = F[2*i];
    a[2*i+1] = F[2*i+1];
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
    if(x[2*i] < XMIN) reflect(XMIN,&x[2*i],&v[2*i],&a[2*i]);
    if(x[2*i] > XMAX) reflect(XMAX,&x[2*i],&v[2*i],&a[2*i]);
    if(x[2*i+1] < YMIN) reflect(YMIN,&x[2*i+1],&v[2*i+1],&a[2*i+1]);
    if(x[2*i+1] > YMAX) reflect(YMAX,&x[2*i+1],&v[2*i+1],&a[2*i+1]);
  }
}


void compute_forces(int n, float* x, float* F)
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
      if(i!=j){
        dx = x[2*j] - x[2*i];
        dy = x[2*j+1] - x[2*i+1];
        lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
        F[2*i] += lj_scalar * dx;     		    // pos account for the direction of the vector from base molecule
        F[2*i+1] += lj_scalar * dy;
        F[2*j] -= lj_scalar * dx;              // neg account for the direction of the vector from non-base molecule
        F[2*j+1] -= lj_scalar * dy;
      }
    }
  }
}

int main(int argc, char** argv)
{
  int npart,i;
  params param;
  param.npart = N_BODY_NUM;
  param.dt = DT;
  param.eps_lj = EPS;
  param.sig_lj = SIG;

  mols mol;

  mol.x = malloc(2*param.npart * sizeof(float));
  mol.v = malloc(2*param.npart * sizeof(float));
  mol.a = malloc(2*param.npart * sizeof(float));
  mol.F = malloc(2*param.npart * sizeof(float));

  npart = init_particles(param.npart, mol.x , mol.v, &param);
  if(npart < param.npart)
  {
    fprintf(stderr, "Could not generate %d particles, Trying %d particles instead\n",param.npart,npart);
    param.npart = npart;
  }
  init_particles_va( param.npart, mol.v,mol.a, &param);
  compute_forces(param.npart,mol.x,mol.F);
  for(i=0;i<ITERS;i++)
  {
    //memset(mol.F, 0 , 2*param.npart * sizeof(float));
    printf("acc10: %f\n", mol.a[i]);
    verletInt1(param.npart,param.dt , mol.x, mol.v,mol.a);
    box_reflect(param.npart,mol.x,mol.v,mol.a );
    compute_forces(param.npart,mol.x,mol.F);
    verletInt2(param.npart,param.dt, mol.x, mol.v, mol.a);
    //computeAcc(param.npart,mol.F,mol.a);
    memset(mol.F, 0 , 2*param.npart * sizeof(float));
    //printf("acc10: %f\n", mol.a[i]);
  }
  for(i=0; i < param.npart; i++)
  {
    printf("nBody-Num: %d Posx: %f Velx: %f Accx: %f Forcex: %f\n",i,
        mol.x[2*i],mol.v[2*i],mol.a[2*i],mol.F[2*i]);
    printf("nBody-Num: %d Posy: %f Vely: %f Accy: %f Forcey: %f\n",i,
        mol.x[2*i+1],mol.v[2*i+1],mol.a[2*i+1],mol.F[2*i+1]);
  }

  free(mol.x);
  free(mol.v);
  free(mol.a);
  free(mol.F);
  printf("Done\n");
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
