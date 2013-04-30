/// stuff happening 
// nvall -o mdCuda mdCuda.cu -g -G 


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define OPTION			2
#define ITERS			10
#define DELTA			128
#define BASE			2048
#define NUM_THREADS_PER_BLOCK 	256
#define NUM_BLOCKS 		1
#define PRINT_TIME 		1
#define TOL			1e-6
#define OMEGA			1.9
#define MAXVAL			10.0
#define MINVAL			-10.0

// cool way to print	    printf("iter = %4d / %-4d \r", i+1, ITERS); fflush(stdout);

// Just a prototype and declaration
//void initializeArray1D(float *arr, int len, int seed);
void initializeArray2D(float *arr, int len, int seed);
void mmm_kij(float* a0, float* b0, float* c0, int length);
struct timespec diff(struct timespec start, struct timespec end);
void cudaErrorCheck(cudaError_t err);

// MMM with global memory

__global__ void kernel_MMM(float* A, float* B, float* C, int len)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;

	float accum = 0;
	int k = 0;
	for(k = 0; k < len; k++)
	{
	    accum += A[i*len+k] * B[(len*k)+j];
	}
	C[i*len+j] = accum;    
}

int main(int argc, char **argv){
    struct timespec time1,time2;
    struct timespec time_stamp;
    // initialization of the necassary vaiables
    //int i,errorNum=0;
    int arrLen = BASE;
    size_t size = arrLen * arrLen * sizeof(float);
    cudaError_t err = cudaSuccess;

    // Timing related variables
    float elapsed_gpu[OPTION];

    // Arrays on GPU global memory
    float *d_A;
    float *d_B;
    float *d_C;

    // Arrays on the host memory
    float *h_A;
    float *h_B;
    float *h_C;
    float *h_base;


    printf("Length of the array = %d bytes\nLength of array = %d (float)\nNumber of iterations = %d\n", size ,arrLen*arrLen, ITERS);

    // Allocate arrays on host memory
    h_A = (float *) malloc(size);
    h_B = (float *) malloc(size);
    h_C = (float *) malloc(size);
    h_base = (float *) malloc(size);

    cudaEvent_t start1,stop1;

    err=cudaEventCreate(&start1);
    cudaErrorCheck(err);

    err = cudaThreadSynchronize();
    cudaErrorCheck(err);
    
    err = cudaEventRecord(start1,0);
    cudaErrorCheck(err);

    printf("About to cudaMalloc\n");
    err = cudaMalloc((void**) &d_A, size);
    cudaErrorCheck(err);

    err = cudaMalloc((void**) &d_B, size);
    cudaErrorCheck(err);

    err = cudaMalloc((void**) &d_C, size);
    cudaErrorCheck(err);
    printf("Finished the cudaMalloc\n");

    // Initialize the host arrays
    printf("\nInitializing the arrays ...");
    initializeArray2D(h_A, arrLen, 2453);
    initializeArray2D(h_B, arrLen, 923874);
    printf("\t... done\n\n");

    // Transfer the arrays to the GPU memory
    err = cudaMemcpy(d_A, h_A, size , cudaMemcpyHostToDevice);
    cudaErrorCheck(err);

    err = cudaMemcpy(d_B, h_B, size , cudaMemcpyHostToDevice);
    cudaErrorCheck(err);

    // create timer and start the timer
    cudaEvent_t start,stop;
    
    err=cudaEventCreate(&start);
    cudaErrorCheck(err);
    
    err = cudaThreadSynchronize();
    cudaErrorCheck(err);
    
    err = cudaEventRecord(start,0);
    cudaErrorCheck(err);



    /// gives the ceiling function for # of blocks  -->  Launch the kernel
    int blocksPerGrid = (arrLen+15)/16;	
    dim3 dimGrid(blocksPerGrid,blocksPerGrid);		
    dim3 dimBlock(16,16);

    // Generate actual cuda call
    printf("Making call to kernel\n");
    kernel_MMM<<< dimGrid,dimBlock >>>(d_A,d_B,d_C,arrLen);
    
    // Check if kernel execution generated an error
    err = cudaGetLastError();
    cudaErrorCheck(err);

    // Transfer the results back to the host
    printf("Waiting for computation to complete...\n");
    
    err = cudaMemcpy( h_C , d_C , size ,cudaMemcpyDeviceToHost);    
    cudaErrorCheck(err);
    
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

    clock_gettime(CLOCK_REALTIME, &time1);
    mmm_kij(h_A,h_B,h_base,arrLen);
    clock_gettime(CLOCK_REALTIME, &time2);
    time_stamp = diff(time1,time2);

    printf("\nFinal times\n");
    printf("ArraySize, GPU time (msec)\n");
    //printf("GPU time: %f (msec)\t Array Size: %d\n", elapsed_gpu[i],BASE+DELTA*i);
    printf("Time to run kernel: %f\n",elapsed_gpu[0]);
    printf("Time to run kernel with memcopy: %f\n",elapsed_gpu[1]);
    printf("Time to run serial code: %lf\n",time_stamp.tv_sec + time_stamp.tv_nsec/1e9);
                                                  


    err = cudaFree(d_A);
    cudaErrorCheck(err);
    
    err = cudaFree(d_B);
    cudaErrorCheck(err);

    err = cudaFree(d_C);
    cudaErrorCheck(err);

    free(h_A);
    free(h_B);
    free(h_C);
    printf("We actually did it! \n");

    return EXIT_SUCCESS;
}

void initializeArray2D(float *arr, int len, int seed) {
    int i;
    float randNum;
    srand(seed);
    int length = len*len;
    for (i = 0; i < length; i++) {
	randNum = (float) rand()/RAND_MAX;
	arr[i] = -10 + randNum * (20);
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

void cudaErrorCheck(cudaError_t err)
{
    if(err!=cudaSuccess)
    {
	fprintf(stderr, "Failed to create timer (error code %s) \n",cudaGetErrorString(err));
	exit(EXIT_FAILURE);
    }
}

