/// stuff happening 
// nvall -o mdCuda mdCuda.cu -g -G -lrt -lm


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#define LINUX       1			    // is this on a linux machine??
#define NEAREST	    0                       // Are we going to use nearest algorithm
#define OPT	    0                       // N^2 or optimized code??

#define EPS	      1
#define SIG	      1e-2
#define CUT	      2.5
#define RCUT	      (CUT*SIG)
#define CUT2	      CUT*CUT
#define PI            3.14159265
#define DT	      0.001           //  0.001 second time increments	     definitely want to change this 
#define N_BODY_NUM    8192
#define XMAX	      (BOX_SIZE/2.0)
#define XMIN	      -(BOX_SIZE/2.0)
#define YMAX	      (BOX_SIZE/2.0)
#define YMIN	      -(BOX_SIZE/2.0)
#define T0	          1
#define MAX_TRIALS    100
#define ITERS         100
#define BOX_SIZE      10.0
#define GRID_NUM      ((BOX_SIZE)/(RCUT))


#define BLOCK_LENGTH(GRID_NUM,BOX_SIZE) (BOX_SIZE/GRID_NUM)		    // size of block that contains GRID_BLOCK_NUM
#define EST_NUM(GRID_NUM,N_BODY_NUM) (N_BODY_NUM/(GRID_NUM*GRID_NUM))

typedef struct sim_param_t {
    int npart;
    float dt;
    float eps_lj;
    float sig_lj;
}params;

typedef struct molecule_t {
    float* x;
    float* v;
    float* a;
    float* F;
}mols; 

__device__ void compute_forces_naive(int n, int k, float* x, float* F);
__device__ void box_reflect(int k, float* x, float* v, float* a);
__device__ void reflect(float wall, float* x, float* v, float* a);
__device__ void verletInt2(int k, float dt, float* x, float* v, float* a);	 
__device__ float compute_LJ_Scalar(float r2, float eps, float sig2);
__device__ void verletInt1(int k, float dt, float* x, float* v, float* a);	    
int init_particles(int n, float* x, float* v, params* param);
void init_particles_va(int n, float* v,float* a, params* param);

// Just a prototype and declaration
struct timespec diff(struct timespec start, struct timespec end);
void cudaErrorCheck(cudaError_t err);

__global__ void kernel_VanDerWaals(float* x, float* v, float* a, float* F, int particles)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    const int iter = ITERS;
    int r,k;
    float dt = 0.0001;
//    for(r=0; r < iter; r++)
//    {
	for(k = 0; k < particles; k++)
	{
	    if(i!=k)
	    {
		verletInt1(k, dt, x, v, a);
		box_reflect(k, x, v, a);
		compute_forces_naive(i, k, x, F);
		verletInt2(k, dt, x, v, a);
	    
	    }
	}
	memset(F,0,2*particles*sizeof(float));
	__syncthreads();
  //  }
}


int main(int argc, char **argv){
    
    int nsize = N_BODY_NUM;
    if(argc==2)
    {
	nsize = atoi(argv[1]);	
    }                         
    struct timespec time1,time2;
    struct timespec time_stamp;
    cudaError_t err = cudaSuccess;

    // start timing of entire program
    clock_gettime(CLOCK_REALTIME, &time1);
    // Timing related variables
    float elapsed_gpu[2];

    // global information
    params param;
    param.npart = nsize;
    param.dt = DT;
    param.eps_lj = EPS;
    param.sig_lj = SIG;

    // declare size in bytes
    size_t size = 2 * (param.npart) * sizeof(float);

    // Arrays on GPU global memory
    mols d_mol;

    // Arrays on the host memory
    mols h_mol;

    // Allocate arrays on host memory
    h_mol.x = (float *) malloc(size);
    h_mol.v = (float *) malloc(size);
    h_mol.a = (float *) malloc(size);
    h_mol.F = (float *) malloc(size);

    cudaEvent_t start1,stop1;

    err=cudaEventCreate(&start1);
    cudaErrorCheck(err);

    err = cudaThreadSynchronize();
    cudaErrorCheck(err);

    err = cudaEventRecord(start1,0);
    cudaErrorCheck(err);

    printf("About to cudaMalloc\n");
    err = cudaMalloc((void**) &d_mol.x, size);
    cudaErrorCheck(err);

    err = cudaMalloc((void**) &d_mol.v, size);
    cudaErrorCheck(err);

    err = cudaMalloc((void**) &d_mol.a, size);
    cudaErrorCheck(err);

    err = cudaMalloc((void**) &d_mol.F, size);
    cudaErrorCheck(err);

    printf("Finished the cudaMalloc\n");

    // Initialize the host arrays
    printf("\nInitializing the Particles ...");

    param.npart = init_particles(param.npart, h_mol.x, h_mol.v, &param);
    init_particles_va(param.npart, h_mol.v, h_mol.a, &param);

    printf("\t... done\n\n");

    // Transfer the arrays to the GPU memory
    printf("About to cudaMemcpy\n");
    err = cudaMemcpy(d_mol.x, h_mol.x, size , cudaMemcpyHostToDevice);
    cudaErrorCheck(err);

    err = cudaMemcpy(d_mol.v, h_mol.v, size , cudaMemcpyHostToDevice);
    cudaErrorCheck(err);

    err = cudaMemcpy(d_mol.a, h_mol.a, size , cudaMemcpyHostToDevice);
    cudaErrorCheck(err);

    err = cudaMemcpy(d_mol.F, h_mol.F, size , cudaMemcpyHostToDevice);
    cudaErrorCheck(err);
    printf("Finished cudaMemcpy\n");

    // create timer and start the timer
    cudaEvent_t start,stop;

    err=cudaEventCreate(&start);
    cudaErrorCheck(err);

    err = cudaThreadSynchronize();
    cudaErrorCheck(err);

    err = cudaEventRecord(start,0);
    cudaErrorCheck(err);



    /// gives the ceiling function for # of blocks  -->  Launch the kernel
    int blocksPerGrid = ((param.npart+255)/256);	
    
    printf("\n%d\n",blocksPerGrid);
    
    dim3 dimGrid(blocksPerGrid);		
    dim3 dimBlock(256);

    // Generate actual cuda call
    printf("Making call to kernel\n");
    kernel_VanDerWaals<<< blocksPerGrid,dimBlock >>>(d_mol.x, d_mol.v, d_mol.a, d_mol.F, param.npart);


    // Transfer the results back to the host
    printf("Waiting for computation to complete...\n");
    
    // just added this line for debugging purposes
    err = cudaThreadSynchronize();
    cudaErrorCheck(err);
    // Check if kernel execution generated an error
    err = cudaGetLastError();
    cudaErrorCheck(err);

    err = cudaMemcpy( h_mol.x , d_mol.x , size ,cudaMemcpyDeviceToHost);    
    cudaErrorCheck(err);
    printf("Memcpy #1 complete ...\n");

    err = cudaMemcpy( h_mol.v , d_mol.v , size ,cudaMemcpyDeviceToHost);    
    cudaErrorCheck(err);
    printf("Memcpy #2 complete ...\n");

    err = cudaMemcpy( h_mol.a , d_mol.a , size ,cudaMemcpyDeviceToHost);    
    cudaErrorCheck(err);
    printf("Memcpy #3 complete ...\n");

    err = cudaMemcpy( h_mol.F , d_mol.F , size ,cudaMemcpyDeviceToHost);    
    cudaErrorCheck(err);
    printf("Memcpy #4 complete ...\n");

    printf("Complete!\n");

    // Stop and destroy the timer
    err = cudaThreadSynchronize();
    cudaErrorCheck(err);

    err = cudaEventCreate(&stop1);
    cudaErrorCheck(err);

    err = cudaEventRecord(stop1,0);
    cudaErrorCheck(err);

    err = cudaEventSynchronize(stop1);
    cudaErrorCheck(err);

    err = cudaEventCreate(&stop);
    cudaErrorCheck(err);

    err = cudaEventRecord(stop,0);
    cudaErrorCheck(err);

    err = cudaEventSynchronize(stop);
    cudaErrorCheck(err);

    // Get the time
    cudaEventElapsedTime(&elapsed_gpu[0],start,stop);           // inlcuding kernel call only
    cudaEventElapsedTime(&elapsed_gpu[1],start1,stop1);		// including memcopy

    // Clean up our mess
    err = cudaEventDestroy(start);
    cudaErrorCheck(err);

    err = cudaEventDestroy(stop);
    cudaErrorCheck(err);

    err = cudaEventDestroy(start1);
    cudaErrorCheck(err);

    err = cudaEventDestroy(stop1);
    cudaErrorCheck(err);

    clock_gettime(CLOCK_REALTIME, &time2);
    time_stamp = diff(time1,time2);

    printf("\nFinal times\n");
    printf("ArraySize, GPU time (msec)\n");
    //printf("GPU time: %f (msec)\t Array Size: %d\n", elapsed_gpu[i],BASE+DELTA*i);
    printf("Time to run kernel: %f\n",elapsed_gpu[0]);
    printf("Time to run kernel with memcopy: %f\n",elapsed_gpu[1]);
    printf("Time to run serial code: %lf\n",time_stamp.tv_sec + time_stamp.tv_nsec/1e9);




    err = cudaFree(d_mol.x);
    cudaErrorCheck(err);
    err = cudaFree(d_mol.v);
    cudaErrorCheck(err);
    err = cudaFree(d_mol.a);
    cudaErrorCheck(err);
    err = cudaFree(d_mol.F);
    cudaErrorCheck(err);




    free(h_mol.x);
    free(h_mol.v);
    free(h_mol.a);
    free(h_mol.F);
    printf("We actually did it! \n");

    return EXIT_SUCCESS;
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

void cudaErrorCheck(cudaError_t err)
{
    if(err!=cudaSuccess)
    {
	fprintf(stderr, "cudaError (error code %s) \n",cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
}


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


__device__ float compute_LJ_Scalar(float r2, float eps, float sig2)
{
    if(r2 < (CUT2 * sig2))	    // 
    {
	float frac2 = sig2/r2;
	float frac6 = frac2*frac2*frac2;
	return 24.0*eps/r2 * frac6 *(1.0-2.0*frac6);		
    }
    return 0;
}

__device__ void verletInt1(int k, float dt, float* x, float* v, float* a)	    
{                            
    int two_i = 2*k;                                      // assumes that we havbe 2D data
    v[two_i] = a[two_i] * (dt/2.0);		    // spltwo_it up for a 2D 
    v[two_i+1] = a[two_i+1] * (dt/2.0);
    x[two_i] = v[two_i] * dt;
    x[two_i+1] = v[two_i+1] * dt;
}

__device__ void verletInt2(int k, float dt, float* x, float* v, float* a)	 
{
    int two_i = 2*k;
    int v0 = v[two_i];
    int v1 = v[two_i+1];
    v[two_i] = a[two_i] * dt/2.0;		    // spltwo_it up for 2D
    v[two_i+1] = a[two_i+1] * dt/2.0;
    a[two_i] += (v[two_i]-v0)/dt;
    a[two_i+1] += (v[two_i+1]-v1)/dt;
}

// should check for reflection inbetween 
__device__ void reflect(float wall, float* x, float* v, float* a)
{
    //printf("reflected!");
    *x = (2*wall-(*x));
    *v = -(*v);
    *a = -(*a);
}

__device__ void box_reflect(int k, float* x, float* v, float* a)
{
    int two_i = 2*k;
    if(x[two_i] < XMIN) reflect(XMIN,&x[two_i],&v[two_i],&a[two_i]);
    if(x[two_i] > XMAX) reflect(XMAX,&x[two_i],&v[two_i],&a[two_i]);
    if(x[two_i+1] < YMIN) reflect(YMIN,&x[two_i+1],&v[two_i+1],&a[two_i+1]);
    if(x[two_i+1] > YMAX) reflect(YMAX,&x[two_i+1],&v[two_i+1],&a[two_i+1]);
}


// now only executes over the bodies declared!
// n is me and k is them
__device__ void compute_forces_naive(int n, int k, float* x, float* F)
{
    float eps = EPS;
    float sig = SIG;
    float sig2 = sig*sig;
    float dx,dy,lj_scalar;

    dx = x[2*k] - x[2*n];
    dy = x[2*k+1] - x[2*n+1];
    lj_scalar = compute_LJ_Scalar(dx*dx+dy*dy,eps,sig2);
    F[2*n] += lj_scalar * dx;     		    // pos account for the direction of the vector from base molecule
    F[2*n+1] += lj_scalar * dy;
} 



