
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>

#include <omp.h>


#include <stdio.h>


using namespace std;

class PointBin
{
public:

	void doBinning(string filename, int binDim, int* binArray, int algNumber)
	{

		//open file input stream to check size
		ifstream binFile;

		binFile.open(filename, ios::binary | ios::ate);

		//Get file size
		int numVals = (binFile.tellg() / (2 * sizeof(float)));

		//set number of threads that will be used
		int numThreads = 16;
		omp_set_num_threads(numThreads);

		//init counter to check number of points read in
		int numPoints = 0;

		binFile.close();

		//Read in values

		//checked
		//Split file (Sequential fragment) paralell
		if (algNumber == 3)
		{

			//parallel version
#pragma omp parallel for num_threads(numThreads) shared(numPoints)
			for (int i = 0; i < numThreads; i++) {

				//Calculating the start points of each thread
				int start = (numVals*i + (numThreads - 1)) / numThreads;
				int end = (numVals*(i + 1) + (numThreads - 1)) / numThreads;

				//open file stream
				ifstream infile;

				infile.open(filename, ios::binary | ios::in);

				//seek starting pos
				infile.seekg(start*sizeof(float));


				for (int j = start; j < end; j++) {


					float x;
					float y;

					//Read in point values
					infile.read((char *)&x, sizeof(float));

					infile.read((char *)&y, sizeof(float));



					//classify which bin to increment

					int xBin = (int)floor(x * binDim);

					int yBin = (int)floor(y * binDim);


					if (xBin == binDim)
						xBin--;
					if (yBin == binDim)
						yBin--;


#pragma omp atomic 
					binArray[xBin * binDim + yBin]++;
#pragma omp atomic 
					numPoints++;


				}
				infile.close();
			}


		}
		else if (algNumber == 0)
		{
			//checked
			//Sequential version
			for (int i = 0; i < numThreads; i++) {

				int start = (numVals*i + (numThreads - 1)) / numThreads;
				int end = (numVals*(i + 1) + (numThreads - 1)) / numThreads;

				//open file
				ifstream infile;

				infile.open(filename, ios::binary | ios::in);

				//seek starting pos
				infile.seekg(start*sizeof(float));

				for (int j = start; j < end; j++) {

					float x;
					float y;

					//read in points
					infile.read((char *)&x, sizeof(float));
					infile.read((char *)&y, sizeof(float));

					//classify which bin to increment
					int xBin = (int)floor(x * binDim);
					int yBin = (int)floor(y * binDim);


					if (xBin == binDim)
						xBin--;
					if (yBin == binDim)
						yBin--;

					binArray[xBin * binDim + yBin]++;

					numPoints++;

				}
				infile.close();
			}

		}
		else if (algNumber == 2)
		{

			//nested parallel

			const int threadsPerGroup = 64;

			omp_set_nested(1);

#pragma omp parallel for num_threads(numThreads) shared(numPoints)
			for (int i = 0; i < numThreads; i++) {

				int start = (numVals*i + (numThreads - 1)) / numThreads;
				int end = (numVals*(i + 1) + (numThreads - 1)) / numThreads;

				//open file
				ifstream infile;

				infile.open(filename, ios::binary | ios::in);

				//seek starting pos
				infile.seekg(start*sizeof(float));

#pragma omp parallel for num_threads(threadsPerGroup) shared(i, start, end, infile) schedule(dynamic,1)
				for (int j = start; j < end; j++) {


					float x;
					float y;

#pragma omp critical
					{
						infile.read((char *)&y, sizeof(float));

						infile.read((char *)&x, sizeof(float));
					}


					//classify which bin to increment

					int xBin = (int)floor(x * binDim);

					int yBin = (int)floor(y * binDim);


					if (xBin == binDim)
						xBin--;
					if (yBin == binDim)
						yBin--;


#pragma omp atomic 
					binArray[xBin * binDim + yBin]++;
#pragma omp atomic 
					numPoints++;




				}
				infile.close();
			}
			
		}
		else if (algNumber == 1)
		{
			//sketch
			//Hail mary parralel

			ifstream infile;

			infile.open(filename, ios::binary | ios::in);

#pragma omp parallel for num_threads(8)
			for (int a = 0; a < numVals; a++)
			{


				float x;
				float y;

#pragma omp critical
				{
					//Read in point values
					infile.read((char *)&x, sizeof(float));

					infile.read((char *)&y, sizeof(float));
				}



				//classify which bin to increment

				int xBin = (int)floor(x * binDim);

				int yBin = (int)floor(y * binDim);


				if (xBin == binDim)
					xBin--;
				if (yBin == binDim)
					yBin--;


#pragma omp atomic 
				binArray[xBin * binDim + yBin]++;
#pragma omp atomic 
				numPoints++;

			}
			infile.close();
		}


		//output number of elements read in
		cout << "Number of Points: " << numPoints << endl;




	}
};

