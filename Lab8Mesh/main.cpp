#include <mpi.h>
#include <chrono>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

void transpose(int*& matrix, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			int tmp = matrix[i * size + j];
			matrix[i * size + j] = matrix[j * size + i];
			matrix[j * size + i] = tmp;
		}
	}
}

void matrixMultiplicationMPI(int*& A, int*& B, int*& C, int size, int procRank, int procAmount)
{

	int dim = size;
	int i, j, k, p, ind;
	int temp;

	MPI_Status Status;

	int procPartsize = dim / procAmount;
	cout << procPartsize << endl;
	int procPartElem = procPartsize * dim;

	int* bufA = new int[procPartElem];
	int* bufB = new int[procPartElem];
	int* bufC = new int[procPartElem];

	if (procRank == 0)
		transpose(B, size);

	MPI_Scatter(A, procPartElem, MPI_INT, bufA, procPartElem, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(B, procPartElem, MPI_INT, bufB, procPartElem, MPI_INT, 0, MPI_COMM_WORLD);

	temp = 0;

	for (i = 0; i < procPartsize; i++)
	{
		for (j = 0; j < procPartsize; j++)
		{
			for (k = 0; k < dim; k++)
				temp += bufA[i * dim + k] * bufB[j * dim + k];

			bufC[i * dim + j + procPartsize * procRank] = temp;
			temp = 0;
		}
	}

	int nextProc;
	int prevProc;

	for (p = 1; p < procAmount; p++)
	{
		nextProc = procRank + 1;

		if (procRank == procAmount - 1)
			nextProc = 0;

		prevProc = procRank - 1;

		if (procRank == 0)
			prevProc = procAmount - 1;

		MPI_Sendrecv_replace(bufB, procPartElem, MPI_INT, nextProc, 0, prevProc, 0, MPI_COMM_WORLD, &Status);

		temp = 0;

		for (i = 0; i < procPartsize; i++)
			for (j = 0; j < procPartsize; j++)
			{

				for (k = 0; k < dim; k++)
					temp += bufA[i * dim + k] * bufB[j * dim + k];

				if (procRank - p >= 0)
					ind = procRank - p;

				else
					ind = procAmount - p + procRank;

				bufC[i * dim + j + ind * procPartsize] = temp;
				temp = 0;

			}

	}

	MPI_Gather(bufC, procPartElem, MPI_INT, C, procPartElem, MPI_INT, 0, MPI_COMM_WORLD);

	delete[]bufA;
	delete[]bufB;
	delete[]bufC;
}

int main()
{
	int rank, procAmount;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &procAmount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	srand(time(0));		//init random

	vector<int> sizes = { 100, 1000, 5000 };

	for (int size : sizes) 
	{
		int* A;
		int* B;
		int* C;

		chrono::time_point<std::chrono::high_resolution_clock> start, end;

		if (rank == 0)
		{
			A = new int[size * size];
			B = new int[size * size];
			C = new int[size * size];

			for (int i = 0; i < size * size; i++) {			//init matrices
				A[i] = rand() % 10 + 1;
				B[i] = rand() % 10 + 1;
			}

		}

		start = chrono::high_resolution_clock::now();	//get start time point of multiplication

		matrixMultiplicationMPI(A, B, C, size, rank, procAmount);

		end = chrono::high_resolution_clock::now();

		int millisec = chrono::duration_cast<chrono::milliseconds>(end - start).count();

		if (rank == 0)
		{
			cout << "Time for size = " << size << " is " << millisec << "ms\n";
		}
	}

	MPI_Finalize();
	return 0;
}