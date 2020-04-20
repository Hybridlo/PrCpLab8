#include <chrono>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

void matrixMultiplication(vector<vector<int>>& A, vector<vector<int>>& B, vector<vector<int>>& C, int size)
{
	for (int i = 0; i < size; i++) {
		C.push_back(vector<int>());

		for (int j = 0; j < size; j++) {
			C[i].push_back(0);

			for (int k = 0; k < size; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}
}

int main()
{
	srand(time(0));		//init random

	vector<int> sizes = { 100, 1000, 5000 };
	
	vector<vector<int>> A, B, C;

	for (int size : sizes) {

		for (int i = 0; i < size; i++) {			//init matrices
			A.push_back(vector<int>());
			B.push_back(vector<int>());

			for (int j = 0; j < size; j++) {
				A[i].push_back(rand() % 10 + 1);	//value from 1 to 10
				B[i].push_back(rand() % 10 + 1);
			}
		}

		chrono::time_point<std::chrono::high_resolution_clock> start, end;

		start = chrono::high_resolution_clock::now();	//get start time point of multiplication

		matrixMultiplication(A, B, C, size);

		end = chrono::high_resolution_clock::now();

		int millisec = chrono::duration_cast<chrono::milliseconds>(end - start).count();

		cout << "Time for size = " << size << " is " << millisec << "ms\n";

		A.clear();
		B.clear();
		C.clear();
	}
}