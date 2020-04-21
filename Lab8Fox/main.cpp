#include <mpi.h>
#include <chrono>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

typedef struct
{
	int procAmount;
	MPI_Comm gridComm;
	MPI_Comm rowComm;
	MPI_Comm columnComm;
	int dim;		//size of 2d grid
	int row;
	int column;
	int rank;
} ProcInfo;

void getInfo(ProcInfo* info)
{
	int dimentions[2];
	int isPeriodic[2];
	int coordinates[2];
	int isGetDim[2];

	MPI_Comm_size(MPI_COMM_WORLD, &(info->procAmount));

	info->dim = (int)sqrt((double)info->procAmount);

	dimentions[0] = dimentions[1] = info->dim;

	isPeriodic[0] = isPeriodic[1] = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dimentions, isPeriodic, 1, &(info->gridComm));

	MPI_Comm_rank(info->gridComm, &(info->rank));

	MPI_Cart_coords(info->gridComm, info->rank, 2, coordinates);

	info->row = coordinates[0];
	info->column = coordinates[1];

	isGetDim[0] = 0;
	isGetDim[1] = 1;	//get row comm

	MPI_Cart_sub(info->gridComm, isGetDim, &(info->rowComm));

	isGetDim[0] = 1;	//get column comm
	isGetDim[1] = 0;

	MPI_Cart_sub(info->gridComm, isGetDim, &(info->columnComm));
}

void reorganizeMatrix(int*& matrix, int*& remadeMatrix, int smallSize, int size)	//reorganize matrix so that each matrix goes one after the other
{
	int remadeMatrixPtr = 0;

	for (int ptr = 0; ptr < size; ptr += smallSize)
	{
		for (int i = 0; i < smallSize; i++)
		{
			for (int j = 0; j < smallSize; j++)
			{
				remadeMatrix[remadeMatrixPtr] = matrix[ptr + i * size + j];
				remadeMatrixPtr++;
			}
		}
	}
}

void testReorganizeMatrix()		//to check if it works
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int* m = new int[18] {1, 2, 3, 10, 11, 12,		//2 matrices 3x3
		4, 5, 6, 13, 14, 15,
		7, 8, 9, 16, 17, 18};

	int* r = new int[18];
	reorganizeMatrix(m, r, 3, 6);
	if (rank == 0)
	{
		for (int i = 0; i < 18; i++)
			cout << r[i] << " ";	//should be an array of numbers 1 - 18
		cout << endl;
	}
	
	delete[] m;
	delete[] r;

	m = new int[16]{ 21, 22, 25, 26, 29, 30, 33, 34,		//4 matrices 2x2
					 23, 24, 27, 28, 31, 32, 35, 36 };

	r = new int[16];

	reorganizeMatrix(m, r, 2, 8);
	if (rank == 0)
	{
		for (int i = 0; i < 16; i++)
			cout << r[i] << " ";	//should be an array of numbers 21 - 36
	}

	delete[] m;
	delete[] r;
}

void scatterMatrix(int*& A, int*& B, int*& localA, int*& localB, int size, ProcInfo* info)
{
	int submatrixSize = size / sqrt(info->procAmount);
	int tmpMatrixChunk = size * submatrixSize;
	int* tmpMatrixA = new int[tmpMatrixChunk];
	int* tmpMatrixB = new int[tmpMatrixChunk];

	if (info->column == 0)
	{
		MPI_Scatter(A, tmpMatrixChunk, MPI_INT, tmpMatrixA, tmpMatrixChunk, MPI_INT, 0, info->columnComm);
		MPI_Scatter(B, tmpMatrixChunk, MPI_INT, tmpMatrixB, tmpMatrixChunk, MPI_INT, 0, info->columnComm);
	}

	int* tmpMatrixRemadeA = new int[tmpMatrixChunk];
	int* tmpMatrixRemadeB = new int[tmpMatrixChunk];

	if (info->column == 0)
	{
		reorganizeMatrix(tmpMatrixA, tmpMatrixRemadeA, submatrixSize, size);
		reorganizeMatrix(tmpMatrixB, tmpMatrixRemadeB, submatrixSize, size);
	}

	delete[] tmpMatrixA;
	delete[] tmpMatrixB;

	MPI_Scatter(tmpMatrixRemadeA, submatrixSize * submatrixSize, MPI_INT,
		localA, submatrixSize * submatrixSize, MPI_INT, 0, info->rowComm);

	MPI_Scatter(tmpMatrixRemadeB, submatrixSize * submatrixSize, MPI_INT,
		localB, submatrixSize * submatrixSize, MPI_INT, 0, info->rowComm);

	delete[] tmpMatrixRemadeA;
	delete[] tmpMatrixRemadeB;
}

void testScatterMatrix()
{
	ProcInfo info;
	getInfo(&info);

	int* A;
	int* B;

	int size = 4;

	chrono::time_point<std::chrono::high_resolution_clock> start, end;

	if (info.rank == 0)
	{
		A = new int[size * size];
		B = new int[size * size];

		for (int i = 0; i < size * size; i++)
		{
			A[i] = rand() % 10 + 1;
			B[i] = rand() % 10 + 1;
		}

		cout << "A: ";

		for (int i = 0; i < size * size; i++)
		{
			cout << A[i] << " ";
		}

		cout << "\nB: ";

		for (int i = 0; i < size * size; i++)
		{
			cout << B[i] << " ";
		}
	}

	int submatrixSize = size / sqrt(info.procAmount);

	int* localA = new int[submatrixSize * submatrixSize];
	int* localB = new int[submatrixSize * submatrixSize];

	scatterMatrix(A, B, localA, localB, size, &info);

	cout << "\nI am " << info.rank << " my items: A: ";

	for (int i = 0; i < submatrixSize * submatrixSize; i++)
	{
		cout << localA[i] << " ";
	}

	cout << " B: ";

	for (int i = 0; i < submatrixSize * submatrixSize; i++)
	{
		cout << localB[i] << " ";
	}
}

void reverseReorganizeMatrix(int*& matrix, int*& remadeMatrix, int smallSize, int size)
{
	int remadeMatrixPtr = 0;

	for (int ptr = 0; ptr < size; ptr += smallSize)
	{
		for (int i = 0; i < smallSize; i++)
		{
			for (int j = 0; j < smallSize; j++)
			{
				matrix[ptr + i * size + j] = remadeMatrix[remadeMatrixPtr];
				remadeMatrixPtr++;
			}
		}
	}
}

void testReverseReorganizeMatrix()		//to check if it works
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int* m = new int[18];

	int* r = new int[18]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};

	reverseReorganizeMatrix(m, r, 3, 6);

	if (rank == 0)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				cout << m[i * 6 + j] << " ";	//should be 2 matrices 3x3
			}
			cout << endl;
		}
	}

	delete[] m;
	delete[] r;

	m = new int[16];

	r = new int[16]{21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};

	reverseReorganizeMatrix(m, r, 2, 8);
	if (rank == 0)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				cout << m[i * 8 + j] << " ";	//should be 4 matrices 2x2
			}
			cout << endl;
		}
	}

	delete[] m;
	delete[] r;
}

void gatherMatrix(int*& C, int*& localC, int size, ProcInfo* info)		//reverse scatter
{
	int submatrixSize = size / sqrt(info->procAmount);
	int tmpMatrixChunk = size * submatrixSize;

	int* tmpMatrixRemadeC = new int[tmpMatrixChunk];

	MPI_Gather(localC, submatrixSize * submatrixSize, MPI_INT,
		tmpMatrixRemadeC, submatrixSize * submatrixSize, MPI_INT, 0, info->rowComm);

	int* tmpMatrixC = new int[tmpMatrixChunk];

	if (info->column == 0)
	{
		reverseReorganizeMatrix(tmpMatrixC, tmpMatrixRemadeC, submatrixSize, size);
	}

	delete[] tmpMatrixRemadeC;

	if (info->column == 0)
	{
		MPI_Gather(tmpMatrixC, tmpMatrixChunk, MPI_INT, C, tmpMatrixChunk, MPI_INT, 0, info->columnComm);
	}

	delete[] tmpMatrixC;
}

void testGatherMatrix()
{
	ProcInfo info;
	getInfo(&info);

	int size = 4;

	chrono::time_point<std::chrono::high_resolution_clock> start, end;

	int submatrixSize = size / sqrt(info.procAmount);

	int* localC = new int[submatrixSize * submatrixSize];

	cout << "\nI am " << info.rank << " my items: C: ";

	for (int i = 0; i < submatrixSize * submatrixSize; i++)
	{
		localC[i] = rand() % 10 + 1;
		cout << localC[i] << " ";
	}

	int* C;

	if (info.rank == 0)
	{
		C = new int[size * size];
	}

	gatherMatrix(C, localC, size, &info);

	cout << "\nFinal C: ";

	if (info.rank == 0)
	{
		for (int i = 0; i < size * size; i++)
		{
			cout << C[i] << " ";
		}
	}
}

void matrixMultiply(int*& A, int*& B, int*& C, int size)
{
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			int temp = 0;

			for (int k = 0; k < size; k++)
				temp += A[i * size + k] * B[k * size + j];
			C[i * size + j] = temp;
		}
	}
}

void setToZero(int*& matrix, int size)
{
	for (int i = 0; i < size * size; i++)
		matrix[i] = 0;
}

void matrixMultiplicationMPI(int*& localA, int*& localB, int*& localC, int size, ProcInfo* info)
{
	int* bufA;
	int bcastRoot;
	int submatrixSize;
	int source;
	int destination;
	int tag = 57;
	MPI_Status status;

	submatrixSize = size / info->dim;
	setToZero(localC, submatrixSize);

	source = (info->row + 1) % info->dim;
	destination = (info->row + info->dim - 1) % info->dim;

	bufA = new int[submatrixSize * submatrixSize];

	for (int step = 0; step < info->dim; step++)
	{
		bcastRoot = (info->row + step) % info->dim;

		if (bcastRoot == info->column)
		{
			MPI_Bcast(localA, 1, MPI_INT, bcastRoot, info->rowComm);
			matrixMultiply(localA, localB, localC, submatrixSize);
		}

		else
		{
			MPI_Bcast(bufA, 1, MPI_INT, bcastRoot, info->rowComm);
			matrixMultiply(bufA, localB, localC, submatrixSize);
		}

		MPI_Send(localB, 1, MPI_INT, destination, tag, info->columnComm);

		MPI_Recv(localB, 1, MPI_INT, source, tag, info->columnComm, &status);
	}
}

int main()
{
	MPI_Init(NULL, NULL);

	ProcInfo info;
	getInfo(&info);

	srand(time(0) + info.rank);		//init random + rank for gather test

	if (sqrt(info.procAmount) - (int)sqrt(info.procAmount) != 0)
	{
		if (info.rank == 0)	cout << "Error in process number, aborting";

		MPI_Finalize();
		return 0;
	}

	//testReorganizeMatrix();
	//testReverseReorganizeMatrix();
	//testScatterMatrix();
	//testGatherMatrix();

	vector<int> sizes = { 100, 1000, 5000 };

	for (int size : sizes)
	{
		int* A;
		int* B;
		int* C;

		chrono::time_point<std::chrono::high_resolution_clock> start, end;

		if (info.rank == 0)
		{
			A = new int[size * size];
			B = new int[size * size];
			C = new int[size * size];

			for (int i = 0; i < size * size; i++)			//init matrices
			{
				A[i] = rand() % 10 + 1;
				B[i] = rand() % 10 + 1;
			}

		}

		int submatrixSize = size / sqrt(info.procAmount);

		int* localA = new int[submatrixSize * submatrixSize];
		int* localB = new int[submatrixSize * submatrixSize];
		int* localC = new int[submatrixSize * submatrixSize];

		start = chrono::high_resolution_clock::now();	//get start time point of multiplication

		scatterMatrix(A, B, localA, localB, size, &info);

		matrixMultiplicationMPI(localA, localB, localC, size, &info);

		gatherMatrix(C, localC, size, &info);

		end = chrono::high_resolution_clock::now();

		int millisec = chrono::duration_cast<chrono::milliseconds>(end - start).count();

		if (info.rank == 0)
		{
			cout << "Time for size = " << size << " is " << millisec << "ms\n";
		}
	}

	MPI_Finalize();
	return 0;
}