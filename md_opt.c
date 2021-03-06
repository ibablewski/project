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

#define LINUX       1			    // is this on a linux machine??
#define NEAREST	    1                       // Are we going to use nearest algorithm
#define OPT	    0                       // N^2 or optimized code??

#define EPS	      1
#define SIG	      1e-2
#define CUT	      3
#define RCUT	      (CUT*SIG)
#define CUT2	      CUT*CUT
#define PI            3.14159265
#define DT	      0.0001           //  0.001 second time increments
#define N_BODY_NUM    1000
#define XMAX	      (BOX_SIZE/2.0)
#define XMIN	      -(BOX_SIZE/2.0)
#define YMAX	      (BOX_SIZE/2.0)
#define YMIN	      -(BOX_SIZE/2.0)
#define T0	          1
#define MAX_TRIALS    10
#define ITERS         100
#define BOX_SIZE      10.0
#define GRID_NUM      ((BOX_SIZE)/(RCUT))

#define BLOCK_LENGTH(GRID_NUM,BOX_SIZE) (BOX_SIZE/GRID_NUM)		    // size of block that contains GRID_BLOCK_NUM
#define EST_NUM(GRID_NUM,N_BODY_NUM) (N_BODY_NUM/(GRID_NUM*GRID_NUM))


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

typedef struct adjacents {
  int* n;
}adjacent;



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
    v[2*i+1] = a[2*i+1] * dt/2.0;
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
  printf("reflected!");
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
 
void compute_forces_naive(int n, float* x, float* F)
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
      //F[2*j] -= lj_scalar * dx;              // neg account for the direction of the vector from non-base molecule
      //F[2*j+1] -= lj_scalar * dy;            // only applies in the non naive serial version
       }
    }
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
    for(j = i+1; j < n; j++)
    {
      // if(i!=j){
      dx = x[2*j] - x[2*i];
      dy = x[2*j+1] - x[2*i+1];
      lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
      F[2*i] += lj_scalar * dx;     		    // pos account for the direction of the vector from base molecule
      F[2*i+1] += lj_scalar * dy;
      F[2*j] -= lj_scalar * dx;              // neg account for the direction of the vector from non-base molecule
      F[2*j+1] -= lj_scalar * dy;            // used in the non naive serial method 
      // }
    }
  }
}

void compute_forces_nearby(int n,int* adj, float* x, float* F, int blockNum1D, int partsPerBlock)
{
  int i,j,k,l,r, two_i;
  float eps = EPS;
  float sig = SIG;
  float sig2 = sig*sig;
  float dx,dy,lj_scalar;
  int myBlock,tmp,index;

  for(i=0;i<n;i++)
  {
    //printf("%d\n",i);
    two_i = 2*i;
    myBlock = getMyBlock(n,two_i,adj,partsPerBlock);
    // work along the top row
    if(myBlock < blockNum1D)
    {
      // work as if on top row
      if(myBlock!=0 && myBlock!=(blockNum1D-1))
      {
        for(r = 0; r <= blockNum1D; r+=blockNum1D)
        {
          for(l=myBlock + r-1 ;l <= myBlock + r +1 ; l++)
          {
            tmp = l * partsPerBlock;
            for(k = 0; k < partsPerBlock; k++)
            {
              if(adj[tmp+k]!=-1)
              {
                index = tmp+k;
                dx = x[index] - x[two_i];
                dy = x[index+1] - x[two_i+1];
                lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
                F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
                F[two_i+1] += lj_scalar * dy;
                k++;
              }  
              else
              {
                break;
              }
            }
          }
        }

        if(myBlock==0)
        {
          // for the left corner
          for(r = 0; r <= blockNum1D; r+=blockNum1D)
          {
            for(l=myBlock + r;l <= myBlock + r +1 ; l++)
            {
              tmp = l * partsPerBlock;
              for(k = 0; k < partsPerBlock; k++)
              {
                if(adj[tmp+k]!=-1)
                {
                  index = tmp+k;
                  dx = x[index] - x[two_i];
                  dy = x[index+1] - x[two_i+1];
                  lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
                  F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
                  F[two_i+1] += lj_scalar * dy;
                  k++;
                }  
                else
                {
                  break;
                }
              }  
            }
          } 
        }

        else if(myBlock==blockNum1D-1)
        {
          // for the right top corner
          for(r = 0; r <= blockNum1D; r+=blockNum1D)
          {
            for(l=myBlock + r-1;l <= myBlock + r; l++)
            {
              tmp = l * partsPerBlock;
              for(k = 0; k < partsPerBlock; k++)
              {
                if(adj[tmp+k]!=-1)
                {
                  index = tmp+k;
                  dx = x[index] - x[two_i];
                  dy = x[index+1] - x[two_i+1];
                  lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
                  F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
                  F[two_i+1] += lj_scalar * dy;
                  k++;
                }  
                else
                {
                  break;
                }
              }  
            }
          } 
        }
      }
    }

    else if(myBlock >= (blockNum1D*(blockNum1D-1)))
    {
      // bottom left corner
      if(myBlock==(blockNum1D*blockNum1D-blockNum1D))
      {
        for(r =(-1)*blockNum1D; r <= 0; r+=blockNum1D)
        {
          for(l=myBlock + r;l <= myBlock + r +1 ; l++)
          {
            tmp = l * partsPerBlock;
            for(k = 0; k < partsPerBlock; k++)
            {
              if(adj[tmp+k]!=-1)
              {
                index = tmp+k;
                dx = x[index] - x[two_i];
                dy = x[index+1] - x[two_i+1];
                lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
                F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
                F[two_i+1] += lj_scalar * dy;
                k++;
              }  
              else
              {
                break;
              }
            }  
          }
        }
      }

      else if(myBlock==(blockNum1D*blockNum1D-1))
      {
        // bottom right corner
        for(r =(-1)*blockNum1D; r <= 0; r+=blockNum1D)
        {
          for(l=myBlock + r-1;l <= myBlock + r; l++)
          {
            tmp = l * partsPerBlock;
            for(k = 0; k < partsPerBlock; k++)
            {
              if(adj[tmp+k]!=-1)
              {
                index = tmp+k;
                dx = x[index] - x[two_i];
                dy = x[index+1] - x[two_i+1];
                lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
                F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
                F[two_i+1] += lj_scalar * dy;
                k++;
              }  
              else
              {
                break;
              }
            }  
          }
        }
      }

      else{
        // work as being on bottom row
        for(r =(-1)*blockNum1D; r <= 0; r+=blockNum1D)
        {
          for(l=myBlock + r-1 ;l <= myBlock + r +1 ; l++)
          {
            tmp = l * partsPerBlock;
            for(k = 0; k < partsPerBlock; k++)
            {
              if(adj[tmp+k]!=-1)
              {
                index = tmp+k;
                dx = x[index] - x[two_i];
                dy = x[index+1] - x[two_i+1];
                lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
                F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
                F[two_i+1] += lj_scalar * dy;
                k++;
              }  
              else
              {
                break;
              }
            }  
          }
        }
      }  
    }

    else if(!(myBlock%blockNum1D)) 
    {
      // work as the leftmost column
      for(r =(-1)*blockNum1D; r <= blockNum1D; r+=blockNum1D)
      {
        for(l=myBlock + r;l <= myBlock + r +1 ; l++)
        {
          tmp = l * partsPerBlock;
          for(k = 0; k < partsPerBlock; k++)
          {
            if(adj[tmp+k]!=-1)
            {
              index = tmp+k;
              dx = x[index] - x[two_i];
              dy = x[index+1] - x[two_i+1];
              lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
              F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
              F[two_i+1] += lj_scalar * dy;
              k++;
            }  
            else
            {
              break;
            }
          }  
        }
      }  
    }

    else if((myBlock%blockNum1D)==(blockNum1D-1)) 
    {
      // work as the rightmost column
      for(r =(-1)*blockNum1D; r <= blockNum1D; r+=blockNum1D)
      {
        for(l=myBlock + r-1 ;l <= myBlock + r +1 ; l++)
        {
          tmp = l * partsPerBlock;
          for(k = 0; k < partsPerBlock; k++)
          {
            if(adj[tmp+k]!=-1)
            {
              index = tmp+k;
              dx = x[index] - x[two_i];
              dy = x[index+1] - x[two_i+1];
              lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
              F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
              F[two_i+1] += lj_scalar * dy;
              k++;
            }  
            else
            {
              break;
            }
          }  
        }
      }  
    }

    else
    {
      // work as a middle block
      for(r =(-1)*blockNum1D; r <= blockNum1D; r+=blockNum1D)
      {
        for(l=myBlock + r-1 ;l <= myBlock + r +1 ; l++)
        {
          tmp = l * partsPerBlock;
          for(k = 0; k < partsPerBlock; k++)
          {
            if(adj[tmp+k]!=-1)
            {
              index = tmp+k;
              dx = x[index] - x[two_i];
              dy = x[index+1] - x[two_i+1];
              lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
              F[two_i] += lj_scalar * dx;     	    // pos account for the direction of the vector from base molecule
              F[two_i+1] += lj_scalar * dy;
              k++;
            }  
            else
            {
              break;
            }
          }
        } 
      } 
    }  
  }
}

int maxNumPartPerBlock(float blockSize, float sig)
{
  return ((blockSize*blockSize)/(sig*sig));
}

void gridSort(int n, int numBlocks, int numPartsPerBox, int* adj, float* x)
{
  float start = (BOX_SIZE/2.0)-BOX_SIZE;
  float step = RCUT;                            // width of each block
  int totBlocks = (int) GRID_NUM*GRID_NUM;       // blocks for both in the x and y directions
  int cnt[totBlocks];                            // total # of blocks
  float xblk, yblk;

  int i,j,k,tmp;

  memset(cnt,0,totBlocks * sizeof(int));

  for(i=0;i<n;i++)                            // for each nbody
  {
    /* Look over the X blocks */
    xblk = start;
    for(j=0; j<numBlocks; j++)
    {
      xblk+=step;
      if(x[2*i] < xblk)
      {
        /* Look over the Y blocks */
        yblk = start;
        for(k=0; k < numBlocks; k++)
        {
          yblk+=step;

          /* Execute only if in this X,Y blocks */
          if(x[2*i+1] < yblk)
          {
            tmp = cnt[j+k*numBlocks]++;
            //printf("j: %d, k: %d, Block: %d, ActualMap: %d\n",j,k, j+k*numBlocks,(j+k*numBlocks)*numPartsPerBox+tmp);
            adj[(j+k*numBlocks)*numPartsPerBox+tmp] = 2*i;

            //printf("(%d,%d) ",(j+k*numB) ,2*i);
            break;
          }
        }
        break;
      }
    }
  }
}

int getMyBlock(int n, int id, int* adj, int numPartsPerBox)
{
  int i = 0;
  int tmp = n*numPartsPerBox;
  for(i=0;i<tmp;i++)
  {
    if(adj[i]==id)
    {
      break;
    }
  }
  return i/numPartsPerBox;
}

int main(int argc, char** argv)
{
  //printf("%f\n",GRID_NUM);
  struct timespec diff(struct timespec start, struct timespec end);
  struct timespec time1, time2;
  struct timespec time_stamp;


  int npart,i,j,Num;

  params param;
  mols mol;
  adjacent adj;

  param.npart = N_BODY_NUM;
  param.dt = DT;
  param.eps_lj = EPS;
  param.sig_lj = SIG;

  mol.x = malloc(2*param.npart * sizeof(float));
  mol.v = malloc(2*param.npart * sizeof(float));
  mol.a = malloc(2*param.npart * sizeof(float));
  mol.F = malloc(2*param.npart * sizeof(float));
#if (NEAREST)
  double Blength = BLOCK_LENGTH(GRID_NUM,BOX_SIZE);
  printf("Blength: %lf\n",Blength);
  Num = EST_NUM(GRID_NUM,N_BODY_NUM);
  Num = 4*Num;
  if(!Num)
  {
     Num = 4;
  }
  
  if(N_BODY_NUM < Num)
  {
    Num = N_BODY_NUM;
  }

  adj.n = malloc(sizeof(int) * GRID_NUM * GRID_NUM * Num);
  memset(adj.n,-1,sizeof(int) * GRID_NUM * GRID_NUM *  Num);
  printf("Num: %d\n",Num);
#endif

  npart = init_particles(param.npart, mol.x , mol.v, &param);
  if(npart < param.npart)
  {
    fprintf(stderr, "Could not generate %d particles, Trying %d particles instead\n",param.npart,npart);
    param.npart = npart;
  }
  else
  {
    fprintf(stdout,"Generated %d particles\n",param.npart);
  }

  init_particles_va( param.npart, mol.v,mol.a, &param);

#if(NEAREST)
  printf("Before gridSort\n");
  gridSort(npart, GRID_NUM, Num, adj.n, mol.x);
  printf("After gridSort\n");
#endif

  /*for(i=0;i<npart;i++)
    printf("myBlockNum: %d\n",getMyBlock(param.npart,2*i, adj.n, Num/2));

    for(i=0; i < param.npart; i++)
    {
    printf("nBody-Num: %d Posx: %f Velx: %f Accx: %f Forcex: %f\n",i,
    mol.x[2*i],mol.v[2*i],mol.a[2*i],mol.F[2*i]);
    printf("nBody-Num: %d Posy: %f Vely: %f Accy: %f Forcey: %f\n",i,
    mol.x[2*i+1],mol.v[2*i+1],mol.a[2*i+1],mol.F[2*i+1]);
    }
    */

#if(LINUX)

  clock_gettime(CLOCK_REALTIME, &time1);

#endif
#if(NEAREST)
  compute_forces_nearby(param.npart, adj.n, mol.x, mol.F, GRID_NUM, Num);

#elif(OPT)  
  compute_forces(param.npart,mol.x,mol.F);
#else
    compute_forces_naive(param.npart,mol.x,mol.F);

#endif
  printf("After First ComputeForces\n");
  for(i=0;i<ITERS;i++)
  {
#if(NEAREST)
    gridSort(npart, GRID_NUM, Num, adj.n, mol.x);		    // Added in the gridSort function for each iteration
#endif
    verletInt1(param.npart,param.dt , mol.x, mol.v,mol.a);
    box_reflect(param.npart,mol.x,mol.v,mol.a );
#if(NEAREST)
    compute_forces_nearby(param.npart, adj.n, mol.x, mol.F, GRID_NUM, Num);
#elif(OPT)
    compute_forces(param.npart,mol.x,mol.F);
#else
    compute_forces_naive(param.npart,mol.x,mol.F);
#endif
    verletInt2(param.npart,param.dt, mol.x, mol.v, mol.a);
    memset(mol.F, 0 , 2*param.npart * sizeof(float));
  }

#if(LINUX)

  clock_gettime(CLOCK_REALTIME, &time2);

#endif

  /*
     for(i=0; i < param.npart; i++)
     {
     printf("nBody-Num: %d Posx: %f Velx: %f Accx: %f Forcex: %f\n",i,
     mol.x[2*i],mol.v[2*i],mol.a[2*i],mol.F[2*i]);
     printf("nBody-Num: %d Posy: %f Vely: %f Accy: %f Forcey: %f\n",i,
     mol.x[2*i+1],mol.v[2*i+1],mol.a[2*i+1],mol.F[2*i+1]);
     }
     */

#if(LINUX)

  double blength = BLOCK_LENGTH(GRID_NUM,BOX_SIZE);
  printf("Boxsize: %lf,Blocksize: %lf,MaxBodiesPerBlock: %d\n",BOX_SIZE,blength, maxNumPartPerBlock(blength,SIG));
  time_stamp = diff(time1,time2);
  printf("Execution time: %lf\n",(double)((time_stamp.tv_sec + (time_stamp.tv_nsec/1.0e9))));

#else

  printf("Boxsize: %lf,BlockNum: %lf,MaxBodiesPerBlock: %d\n",BOX_SIZE,GRID_NUM, maxNumPartPerBlock(GRID_NUM,SIG));

#endif

  /*        // test statements for the grid allocation
            int count =0;
            for(i=0;i<GRID_NUM*GRID_NUM*Num/2;i++)
            {
            if(adj.n[i]!=-1)
            {
            count++;
            }
            }
            printf("%d\n",count);

*/
#if(NEAREST)
  free(adj.n);
#endif
  free(mol.x);
  free(mol.v);
  free(mol.a);
  free(mol.F);
  printf("Done\n");

  return 0;
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


