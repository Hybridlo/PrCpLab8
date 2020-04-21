#include <chrono>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

void matrixMultiplication(int*& A, int*& B, int*& C, int size)
{
	int temp;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			temp = 0;

			for (int k = 0; k < size; k++)
				temp += A[i * size + k] * B[k * size + j];
			C[i * size + j] += temp;
		}
	}
}

int main()
{
	srand(time(0));		//init random

	vector<int> sizes = { 100, 1000, 5000 };
	
	int* A;
	int* B;
	int* C;

	for (int size : sizes) {

		A = new int[size * size];
		B = new int[size * size];
		C = new int[size * size];

		for (int i = 0; i < size; i++) {			//init matrices

			for (int j = 0; j < size; j++) {
				A[i * size + j] = rand() % 10 + 1;	//value from 1 to 10
				B[i * size + j] = rand() % 10 + 1;
				C[i * size + j] = 0;
			}
		}

		chrono::time_point<std::chrono::high_resolution_clock> start, end;

		start = chrono::high_resolution_clock::now();	//get start time point of multiplication

		matrixMultiplication(A, B, C, size);

		end = chrono::high_resolution_clock::now();

		int millisec = chrono::duration_cast<chrono::milliseconds>(end - start).count();

		cout << "Time for size = " << size << " is " << millisec << "ms\n";
	}
}