
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>

//#include <omp.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>;
#include <stdio.h>
#include <device_functions.h>


#include "PointBin.cpp"

using namespace std;

#define BLOCKSIZE 8

//tested
//First attempt - Basic implementation
__global__ void filterGPU(int *bins, int *binsOut, int binDim, int filterSize)
{
	// The x and y index of this thread in the array C
	int cx, cy;

	
	cx = blockIdx.x * blockDim.x + threadIdx.x;

	cy = blockIdx.y * blockDim.y + threadIdx.y;

	//calculate the new value and output it to binsOut
	


	if ((cx < binDim) && (cy < binDim))
	{ 

		__shared__ int overlap;

		overlap = ceil(filterSize / 2.0) -1;

		int* window = (int*)malloc(filterSize*filterSize*sizeof(int));
	
	//copy output back from the 


	int windowCount = 0;

	
			//cycle through all elements

			//gather values within filter window

			windowCount = 0;

			for (int n = -overlap ; n <= overlap; n++)
			{
				for (int m = -overlap ; m <= overlap; m++)
				{
					//check if index is within array boundry

					if (((cx + n) < binDim) && ((cx + n) >= 0) && ((cy + m) < binDim) && ((cy + m) >= 0))
					{
						//add value to window
						window[windowCount] = bins[(cx + n)*binDim + (cy + m)];

						windowCount++;
					}
				}
			}

			//find median 

			//sort window values

			 int temp;

			 for (int k = 0; k < windowCount; k++)
			{
				for (int b = 0; b < windowCount - 1; b++)
				{
					if (window[b + 1] < window[b])
					{
						temp = window[b];

						window[b] = window[b + 1];

						window[b + 1] = temp;
					}
				}
			}

			binsOut[cx*binDim + cy] = window[(int)floor(windowCount / 2.0)];
		
			free(window);

	}
}


//check
//Second attempt - Shared memory for block of input array and overlap (for worst case filterSize 21) (first thread of each block copies over data) (valid for all odd filter size up to 21)
__global__ void filterGPU_SM(int *bins, int *binsOut, int binDim, int filterSize)
{
	// The x and y index of this thread in the array C
	int cx, cy;
	int tx, ty;

	//calculating thread identification values
	cx = blockIdx.x * blockDim.x + threadIdx.x;

	cy = blockIdx.y * blockDim.y + threadIdx.y;

	tx = threadIdx.x;

	ty = threadIdx.y;


	__shared__ unsigned short smem[BLOCKSIZE + 20][BLOCKSIZE + 20];

	//calculate the new value and output it to binsOut

	if ((cx < binDim) && (cy < binDim))
	{
		
		//copy into shared memory
		__shared__ int overlap ;

		overlap = ceil(filterSize / 2.0) - 1;
		//if this is the first thread in the block
		if ((threadIdx.x == 0) && (threadIdx.y == 0))
		{
			//Copy sub array into shared mem
			 
			for (int i = 0; i < (BLOCKSIZE + 2 * (overlap)); i++)
			{
				for (int j = 0; j < (BLOCKSIZE + 2 * (overlap)); j++)
				{
					//copy over value

					if (((cx + (i - overlap) >= 0)) && ((cx + (i - overlap) < binDim)) && ((cy + (j - overlap)) >= 0) && ((cy + (j - overlap)) < binDim))
					{
						smem[i][j] = bins[(cx + (i - overlap))*binDim + (cy + (j - overlap))];
					}
					
				}
			}
		}

		//dynamic array for sorting list
		int* window = (int*)malloc(filterSize*filterSize*sizeof(int));
		


		int windowCount = 0;

		//sync threads to ensure that all values are present 
		__syncthreads();


		//cycle through all elements

		//gather values within filter window

	

		for (int n = -overlap; n <= overlap; n++)
		{
			for (int m = -overlap; m <= overlap; m++)
			{
				//check if index is within array boundry

				if (((cx + n) < binDim) && ((cx + n) >= 0) && ((cy + m) < binDim) && ((cy + m) >= 0))
				{
					//add value to window
					window[windowCount] = smem[tx + n + 10][ty + m + 10]; //bins[(cx + n)*binDim + (cy + m)];

					windowCount++;
				}
			}
		}

		//find median 

		//sort window values

		 int temp;

		for (int k = 0; k < windowCount; k++)
		{
			for (int b = 0; b < windowCount - 1; b++)
			{
				if (window[b + 1] < window[b])
				{
					temp = window[b];

					window[b] = window[b + 1];

					window[b + 1] = temp;
				}
			}
		}

		binsOut[cx*binDim + cy] = window[(int)floor(windowCount / 2.0)];


		free(window);

	}
}


//Attempt 3 - Shared memory for block of input array and overlap handled as global call (each thread copies over its correspoinding value) (valid for all odd filter size up to 21)
__global__ void filterGPU_SM2(int *bins, int *binsOut, int binDim, int filterSize)
{
	// The x and y index of this thread in the array C
	int cx, cy;
	int tx, ty;

	cx = blockIdx.x * blockDim.x + threadIdx.x;

	cy = blockIdx.y * blockDim.y + threadIdx.y;

	tx = threadIdx.x;

	ty = threadIdx.y;


	__shared__ unsigned short smem[BLOCKSIZE][BLOCKSIZE];

	//calculate the new value and output it to binsOut

	if ((cx < binDim) && (cy < binDim))
	{

		//copy into shared memory
		__shared__ int overlap;
		//if this is the first thread in the block
		
			//Copy sub array into shared mem
			overlap = ceil(filterSize / 2.0) - 1;

		smem[tx][ty] = bins[cx*binDim + (cy)];


		__syncthreads();


		int* window = (int*)malloc(filterSize*filterSize*sizeof(int));


		int windowCount = 0;


		//cycle through all elements

		//gather values within filter window

		windowCount = 0;

		for (int n = -overlap ; n <= overlap; n++)
		{
			for (int m = -overlap ; m <= overlap; m++)
			{
				//check if index is within array boundry

				if (((tx + n) < BLOCKSIZE) && ((tx + n) >= 0) && ((ty + m) < BLOCKSIZE) && ((ty + m) >= 0))
				{
					//inside of shared mem

					window[windowCount] = smem[tx + n][ty + m]; //bins[(cx + n)*binDim + (cy + m)];

					windowCount++;
				}
				else
				{
					//outside shared
					if (((cx + n) < binDim) && ((cx + n) >= 0) && ((cy + m) < binDim) && ((cy + m) >= 0))
					{
						//in input so pull from global

						window[windowCount] = bins[(cx + n)*binDim + (cy + m)];

						windowCount++;
					}
				}

				
			}
		}

		//find median 

		//sort window values

		int temp;

		for (int k = 0; k < windowCount; k++)
		{
			for (int b = 0; b < windowCount - 1; b++)
			{
				if (window[b + 1] < window[b])
				{
					temp = window[b];

					window[b] = window[b + 1];

					window[b + 1] = temp;
				}
			}
		}

		binsOut[cx*binDim + cy] = window[(int)floor(windowCount / 2.0)];


		free(window);

	}
}


//Attempt 4 - Shared memory for block of input array and overlap handled as global call, large shared memory sorting space per block (each thread copies over its correspoinding value) ()
__global__ void filterGPU_SM_5(int *bins, int *binsOut, int binDim, int filterSize)
{
	// The x and y index of this thread in the array C
	int cx, cy;
	int tx, ty;

	cx = blockIdx.x * blockDim.x + threadIdx.x;

	cy = blockIdx.y * blockDim.y + threadIdx.y;

	tx = threadIdx.x;

	ty = threadIdx.y;



	//help with sorting
	__shared__ unsigned short sortMem[BLOCKSIZE*BLOCKSIZE][1 * sizeof(short)];

	//help with filtering
	__shared__ unsigned short smem[BLOCKSIZE][BLOCKSIZE];
	

	//calculate the new value and output it to binsOut

	if ((cx < binDim) && (cy < binDim))
	{

		//copy into shared memory
		int overlap = ceil(filterSize / 2.0) - 1;;

		//copy over its value to sm
		smem[tx][ty] = bins[cx*binDim + (cy)];

		int windowCount = 0;

		__syncthreads();


		//cycle through all elements

		//gather values within filter window

		for (int n = -overlap ; n < overlap; n++)
		{
			for (int m = -overlap; m < overlap; m++)
			{
				//check if index is within array boundry

				if (((tx + n) < BLOCKSIZE) && ((tx + n) >= 0) && ((ty + m) < BLOCKSIZE) && ((ty + m) >= 0))
				{
					//inside of shared mem

					sortMem[tx*BLOCKSIZE + ty][windowCount] = smem[tx + n][ty + m]; //bins[(cx + n)*binDim + (cy + m)];

					windowCount++;
				}
				else
				{
					//outside shared
					if (((cx + n) < binDim) && ((cx + n) >= 0) && ((cy + m) < binDim) && ((cy + m) >= 0))
					{
						//in array - pull from global

						sortMem[tx*BLOCKSIZE + ty][windowCount] = bins[(cx + n)*binDim + (cy + m)];

						windowCount++;
					}
				}
			}
		}

		//find median 

		//sort window values

		int temp;

		for (int k = 0; k < windowCount; k++)
		{
			for (int b = 0; b < windowCount - 1; b++)
			{
				if (sortMem[tx*blockDim.x + ty][b + 1] < sortMem[tx*blockDim.x + ty][b])
				{
					temp = sortMem[tx*blockDim.x + ty][b];

					sortMem[tx*blockDim.x + ty][b] = sortMem[tx*blockDim.x + ty][b + 1];

					sortMem[tx*blockDim.x + ty][b + 1] = temp;
				}
			}
		}

		//copy over output value to output array
		binsOut[cx*binDim + cy] = sortMem[tx*blockDim.x + ty][(int)floor(windowCount / 2.0)];

	}
}



//Attempt 4 - Shared memory for block of input array and overlap handled as global call, large shared memory sorting space per block (each thread copies over its correspoinding value) ()
__global__ void filterGPU_SM_F_21(int *bins, int *binsOut, int binDim, int filterSize)
{
	// The x and y index of this thread in the array C
	int cx, cy;
	int tx, ty;

	cx = blockIdx.x * blockDim.x + threadIdx.x;

	cy = blockIdx.y * blockDim.y + threadIdx.y;

	tx = threadIdx.x;

	ty = threadIdx.y;


	//help with filtering
	__shared__ unsigned short smem[BLOCKSIZE + 20][BLOCKSIZE + 20];


	//calculate the new value and output it to binsOut

	if ((cx < binDim) && (cy < binDim))
	{

		//copy into shared memory
		int overlap = ceil(filterSize / 2.0) - 1;

		//first thred load method due to large number of value calls, smem is paramount
		if ((threadIdx.x == 0) && (threadIdx.y == 0))
		{
			//Copy sub array into shared mem

			for (int i = 0; i < (BLOCKSIZE + 2 * (overlap)); i++)
			{
				for (int j = 0; j < (BLOCKSIZE + 2 * (overlap)); j++)
				{
					//copy over value

					if (((cx + (i - overlap) >= 0)) && ((cx + (i - overlap) < binDim)) && ((cy + (j - overlap)) >= 0) && ((cy + (j - overlap)) < binDim))
					{
						smem[i][j] = bins[(cx + (i - overlap))*binDim + (cy + (j - overlap))];
					}

				}
			}
		}

		int windowCount = 0;

		__syncthreads();

		int smallerThan;
		int equalTo;

		int currentPossibleMedian;

		int medianPos = (int)floor(windowCount / 2.0);


		//cycle through all elements

		
		for (int i = -overlap; i < overlap; i++)
		{
			for (int j = -overlap; j < overlap; j++)
			{
				//get current possible mean
				currentPossibleMedian = (int)smem[i*BLOCKSIZE + j];
		

		//iterating over to count number that are smaller and equal

		for (int n = -overlap; n < overlap; n++)
		{
			for (int m = -overlap; m < overlap; m++)
			{
				//check if index is within array boundry
					if (((cx + n) < binDim) && ((cx + n) >= 0) && ((cy + m) < binDim) && ((cy + m) >= 0))
					{
						//in array - pull from global
						//add to group count
						windowCount++;

						if (currentPossibleMedian > smem[n][m])
						{
							//smaller than possible mean
							smallerThan++;
						}

						if (currentPossibleMedian == smem[n][m])
						{
							//equal to the possible mean
							equalTo++;
						}
					}				
			}
		}

		if ((smallerThan < medianPos + 1) && ((smallerThan + equalTo) >= medianPos))
		{
			binsOut[cx*binDim + cy] = smem[i][j];
		}

			}

	}

	
	}
}



int main(int argc, char * argv[])
{


	clock_t startT, stopT;

	cout << "HPC Assignment" << endl;
	cout << "Louis Buys" << endl;

	 const int filtSize = atoi(argv[2]);

	cout << "======================================" << endl << endl;

	

	//Somewhere here I must read in the points and "bin" them apparently


	//args are filename bins 
	if ((argc < 3) != true)
	{
		//has all input

		//string filename = "C:\\Users\\Louis\\Documents\\Visual Studio 2013\\Projects\\HPC_GPU\\Debug\\Points_[1.0e+08]_Noise_[030]_Normal.bin";
		cout << "attempting to access specified file" << endl;

		string filename = argv[1]; // "D:\\Points_[1.0e+08]_Noise_[030]_Normal.bin";

		

		ifstream binFile;

		binFile.open(filename, ios::binary | ios::ate);

		//check if file exists

		if (binFile.is_open() != false)
		{
			cout << "File found" << endl;
		}
		else
		{
			cout << "Missing file." << endl;
			system("pause");
			exit(0);
		}
		

		//check size of file (Number of points to read in)
		int numVals = (binFile.tellg() / (2 * sizeof(float)));


		// Time to perform various tasks [ms]
		float readFilePar, readFileSeq, filterCPUcalTime, filterGPUcalTime;

		//collecting command line input
		int binDim = atoi(argv[3]);
		int filterSize = atoi(argv[2]);

		cout << "Using filter size: " << filterSize << endl;
		cout << "Using binDim size: " << binDim << endl;

		//allocate memory for bin array SEQ
		int* binArrayPAR = NULL;
		int* binArrayOutPAR = NULL;

		//Allocate memory for Bins SEQ
		binArrayPAR = (int*)malloc(binDim*binDim*sizeof(int));
		binArrayOutPAR = (int*)malloc(binDim*binDim*sizeof(int));

		//allocate memory for bin array SEQ
		int* binArraySEQ = NULL;
		int* binArrayOutSEQ = NULL;

		//Allocate memory for Bins SEQ
		binArraySEQ = (int*)malloc(binDim*binDim*sizeof(int));
		binArrayOutSEQ = (int*)malloc(binDim*binDim*sizeof(int));

		//Initialize both arrays to 0
		for (int i = 0; i < binDim; i++)
		{
			for (int j = 0; j < binDim; j++)
			{
				//init all elements to 0
				binArraySEQ[i*binDim + j] = 0;
				binArrayOutSEQ[i*binDim + j] = 0;
				binArrayPAR[i*binDim + j] = 0;
				binArrayOutPAR[i*binDim + j] = 0;

			}
		}

		//final param for dobinning to slect algorithms
		//0 seq
		//1 hail mary
		//2 nested
		//3 split

		//init object
		PointBin *pb = new PointBin();

		//Reading in values from file

		//start timing
		startT = clock();


		pb->doBinning(filename, binDim, binArrayPAR, 2);

		//stop timing
		stopT = clock();

		readFilePar = (stopT - startT)/1000;

		cout << "Read File in Paralell Time: " << readFilePar << " s" << endl;

		//Now read file in seq

		startT = clock();

		pb->doBinning(filename, binDim, binArraySEQ, 0);

		stopT = clock();

		readFileSeq = (stopT - startT)/1000;

		cout << "Read File Sequentially Time: " << readFileSeq << " s" << endl;

		//////////////////////////////////////////////////////////////////
		//Speedup calc

		float readSpeedUp = readFileSeq / readFilePar;

		cout << "File read speedUp = " << readSpeedUp << endl;

		//perform the 2D median filter
		cout << "Performing filter with size : " << filterSize << endl;


		//////////////////////////////////////////////////////////////////////
		//SEQ vs PAR parity checking

		bool eq = true;

		for (int i = 0; i < binDim; i++)
		{
			for (int j = 0; j < binDim; j++)
			{
				if (binArraySEQ[i*binDim + j] != binArrayPAR[i*binDim + j])
				{
					eq = false;
			//		cout << binArraySEQ[i*binDim + j] << "!=" << binArrayPAR[i*binDim + j] << endl;
				}
			}
		}

		if (eq != true)
		{
			cout << "ERROR : PAR and SEQ don't match" << endl;
			//system("pause");
		}
		else
		{
			cout << "SUCCESS : PAR and SEQ match" << endl;
		}
		
		////////////////////////////////////////////////////////////////////

		cout << "============================================" << endl;
		cout << "End of Binning" << endl;
		cout << "============================================" << endl << endl;

		cout << "Starting filtering..." << endl;

		//GPU version

		cout << "GPU version" << endl;

		//start timer
		startT = clock();

			//select device
			cudaSetDevice(0);
	
			//initialize error checking and logging scaffolding
			cudaEvent_t start, end;
			cudaEvent_t memHtoD, memDtoH;
			
			cudaEventCreate(&start);
			cudaEventCreate(&end);
					
			cudaEventCreate(&memHtoD);
			cudaEventCreate(&memDtoH);

			cudaError_t result;

		// Device memory pointers
		int* d_ArrayBins = NULL;         // Device array A
		int* d_ArrayBinsOut = NULL;         // Device array B

		//allocate device memory
		cudaMalloc(&d_ArrayBins, binDim*binDim*sizeof(int));
		cudaMalloc(&d_ArrayBinsOut, binDim*binDim*sizeof(int));

		//record event
		cudaEventRecord(memHtoD);

		//copy data to device
		cudaMemcpy(d_ArrayBins, binArrayPAR, binDim*binDim*sizeof(int), cudaMemcpyHostToDevice);


		result = cudaDeviceSynchronize();   // Wait for the device to finish copying memory
		cudaEventSynchronize(memHtoD);

		if (result != cudaSuccess)
			fprintf(stderr, "\nERROR: Error copying data to device [%s]\n", cudaGetErrorString(result));
		else
			printf("Data copy successful\n");
		

		//////////////////////////////////////////////////////////////////////////////////////

		// Variables to describe the blocks and grid


		dim3 threadsPerBlock ;

		threadsPerBlock.x = BLOCKSIZE;
		threadsPerBlock.y = BLOCKSIZE;

		int numBlocksPerRow = ceil(sqrt((binDim*binDim) / (BLOCKSIZE*BLOCKSIZE)));

		dim3 blocksPerGrid(numBlocksPerRow, numBlocksPerRow, 1);

		cudaEventRecord(start);


		//launch kernel
		filterGPU_SM_F_21 << < blocksPerGrid, threadsPerBlock >> >(d_ArrayBins, d_ArrayBinsOut, binDim, filterSize);



				result = cudaGetLastError();  // This determines whether the kernel was launched

				if (result == cudaSuccess)
					printf("Running kernel\n");
				else
				{
					fprintf(stderr, "ERROR: Error at kernel launch %s\n", cudaGetErrorString(result));
				}
			
				cudaEventRecord(end);
				result = cudaDeviceSynchronize();  // This will return when the kernel computation is complete, remember asynchronous execution

				// Complete message
				if (result == cudaSuccess)
				{
				
					cudaEventSynchronize(end);
			
				}
				else
				{
					system("pause");
					fprintf(stderr, "\nERROR: Error after kernel launch %s\n", cudaGetErrorString(result));
					
				}
					
		//copy data back from device
				cudaEventRecord(memDtoH);


		cudaMemcpy(binArrayOutPAR, d_ArrayBinsOut, binDim*binDim*sizeof(int), cudaMemcpyDeviceToHost);

		result = cudaDeviceSynchronize();   // Wait for the device to finish copying memory
		cudaEventSynchronize(memDtoH);

		if (result != cudaSuccess)
			fprintf(stderr, "\nERROR: Error copying data from device [%s]\n", cudaGetErrorString(result));
		else
			printf("Data copy successful\n");
		
		//free device memory
		cudaFree(d_ArrayBins);
		cudaFree(d_ArrayBinsOut);

		//stop clock
		stopT = clock();

		filterGPUcalTime = (stopT - startT)/1000;

		cout << "Filter GPU Time: " << filterGPUcalTime << " s" << endl;


		//sequential version

		startT = clock();

		//allocate space for median calculation window
		int* window = NULL;
		window = (int*)malloc(filterSize*filterSize*sizeof(int));


		int windowCount = 0;

	
		for (int i = 0; i < binDim; i++)
		{
			for (int j = 0; j < binDim; j++)
			{
				//cycle through all elements

				//gather values within filter window

				windowCount = 0;

				for (int n = -ceil(filterSize / 2.0) + 1; n < ceil(filterSize / 2.0); n++)
				{
					for (int m = -ceil(filterSize / 2.0) + 1; m < ceil(filterSize / 2.0); m++)
					{
						//check if index is within array boundry

						if (((i + n) < binDim) && ((i + n) >= 0) && ((j + m) < binDim) && ((j + m) >= 0))
						{
							//add value to window
							window[windowCount] = binArraySEQ[(i + n)*binDim + (j + m)];

							windowCount++;
						}
					
					}
				}

				//find median 

				//sort window values

				int temp;

				for (int k = 0; k < windowCount; k++)
				{
					for (int b = 0; b < windowCount - 1; b++)
					{
						if (window[b + 1] < window[b])
						{
							temp = window[b];

							window[b] = window[b + 1];

							window[b + 1] = temp;
						}
					}
				}
				//output median value to out array
				binArrayOutSEQ[i*binDim + j] = window[(int)floor(windowCount / 2.0)];

			}
		}

	//stop clock
		stopT = clock();

		filterCPUcalTime = (stopT - startT)/1000;

		cout << "Filter CPU Time: " << filterCPUcalTime << " s" << endl;

		//Clac speedup
		float filterSpeedUp = filterCPUcalTime / filterGPUcalTime;

		cout << "filter speedUp = " << filterSpeedUp << endl;

		//////////////////////////////////////////////////////////////////////
		//SEQ vs PAR parity checking

		 eq = true;

		for (int i = 0; i < binDim; i++)
		{
			for (int j = 0; j < binDim; j++)
			{
				if (binArrayOutSEQ[i*binDim + j] != binArrayOutPAR[i*binDim + j])
				{
					eq = false;
				//	cout << binArraySEQ[i*binDim + j] << "!=" << binArrayPAR[i*binDim + j] << endl;
				}
			}
		}

		if (eq != true)
		{
			cout << "ERROR : PAR and SEQ don't match" << endl;
		//	system("pause");
		}
		else
		{
			cout << "SUCCESS : PAR and SEQ match" << endl;
		}

		////////////////////////////////////////////////////////////////////



		//output csv file

		ofstream myfile("out.csv");
		if (myfile.is_open())
		{


			for (int i = 0; i < binDim; i++)
			{
				myfile << "," << (((i + 1)* (1.0 / binDim)) + (((1.0 / binDim) / 2)));
			}

			myfile << endl;

			for (int i = 0; i < binDim; i++)
			{

				myfile << (((i + 1)* (1.0 / binDim)) + ((1.0 / binDim) / 2));

				for (int k = 0; k < binDim; k++)
				{
					myfile << "," << binArrayOutPAR[i*binDim + k];
				}
				myfile << endl;
			}


			myfile.close();
		}

		//release memory
		binArraySEQ = NULL;
		binArrayOutSEQ = NULL;

		binArrayPAR = NULL;
		binArrayOutPAR = NULL;

		//Graphing output

		cout << numVals << "," << filterSize << "," << binDim << "," << readSpeedUp << "," << filterSpeedUp << endl;


	}
	else
	{
		//missing args
		cout << "Please check arguments are : filename and filer size" << endl;
	}

	system("pause");


	return 0;
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     