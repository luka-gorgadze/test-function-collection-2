#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "cg_user.h"
#include "cg_descent.h"
#include "cg_blas.h"  
#include<chrono>
#include<thread>
#include<algorithm>
#include<numeric>
#include<vector>
#include<iostream>
#include<future>
using namespace std;

#include "testFunction.h"
#include "genericFunctions.h"
#include "testFunctions.h"
#include "parallelTestFunctions.h"

const int REPETITIONS = 10;
const int TESTS_LENGTH = 1;
const int TESTS[] = { 37 };
const int RANGE_START = 100000;
const int RANGE_END = 100000;

/*
ARWHEAD 0 5-7
BDQRTIC 1 5-7
COSINE 5 5-7
DIXON3DQ 22 6-8
DQDRTIC 23 5-7
DQRTIC 24 5-7
EDENSCH 25 5-7
EG2 26 bug?
ENGVAL1 27
FLETCHR 29 5-7
FREUROTH 30 5-7
GENHUMPS 31 5-7
GENROSE 32 6-8
LIARWDH 33 6-8
NONDQUAR 37 4-6
POWER 40 6-8
SROSENBR 42 6-8
*/

void runPerformanceTests() {
	double *x;
	double *g;
	sequential::makeTestCollection(nullptr);
	parallel::makeParallelTestCollection(nullptr);

	for (int index = 0; index < TESTS_LENGTH; index++) {
		cout << sequential::test[TESTS[index]].name << endl;
		for (int i = RANGE_START; i <= RANGE_END; i *= 10) {
			double* sequentialValues = (double *)malloc(REPETITIONS * sizeof(double));
			double* sequentialGrads = (double *)malloc(REPETITIONS * sizeof(double));
			double* sequentialValgrads = (double *)malloc(REPETITIONS * sizeof(double));
			double* parallelValues = (double *)malloc(REPETITIONS * sizeof(double));
			double* parallelGrads = (double *)malloc(REPETITIONS * sizeof(double));
			double* parallelValgrads = (double *)malloc(REPETITIONS * sizeof(double));

			for (int r = 0; r < REPETITIONS; r++) {
				//sequential value
				{
					x = (double *)malloc(i * sizeof(double));
					g = (double *)malloc(i * sizeof(double));
					sequential::test[TESTS[index]].StartingGess(x, i);
					auto st = chrono::high_resolution_clock::now();
					sequential::test[TESTS[index]].value(x, i);
					auto diff = chrono::high_resolution_clock::now() - st;
					auto time = chrono::duration_cast<chrono::milliseconds>(diff);
					sequentialValues[r] = time.count();
					free(x);
					free(g);
				}
				//parallel value
				{
					x = (double *)malloc(i * sizeof(double));
					g = (double *)malloc(i * sizeof(double));
					parallel::paralleltest[TESTS[index]].StartingGess(x, i);
					auto st = chrono::high_resolution_clock::now();
					parallel::paralleltest[TESTS[index]].value(x, i);
					auto diff = chrono::high_resolution_clock::now() - st;
					auto time = chrono::duration_cast<chrono::milliseconds>(diff);
					parallelValues[r] = time.count();
					free(x);
					free(g);
				}

				//sequential grad
				{
					x = (double *)malloc(i * sizeof(double));
					g = (double *)malloc(i * sizeof(double));
					sequential::test[TESTS[index]].StartingGess(x, i);
					auto st = chrono::high_resolution_clock::now();
					sequential::test[TESTS[index]].grad(x, g, i);
					auto diff = chrono::high_resolution_clock::now() - st;
					auto time = chrono::duration_cast<chrono::milliseconds>(diff);
					sequentialGrads[r] = time.count();
					free(x);
					free(g);
				}
				//parallel grad
				{
					x = (double *)malloc(i * sizeof(double));
					g = (double *)malloc(i * sizeof(double));
					parallel::paralleltest[TESTS[index]].StartingGess(x, i);
					auto st = chrono::high_resolution_clock::now();
					parallel::paralleltest[TESTS[index]].grad(x, g, i);
					auto diff = chrono::high_resolution_clock::now() - st;
					auto time = chrono::duration_cast<chrono::milliseconds>(diff);
					parallelGrads[r] = time.count();
					free(x);
					free(g);
				}

				//sequential valgrad
				{
					x = (double *)malloc(i * sizeof(double));
					g = (double *)malloc(i * sizeof(double));
					sequential::test[TESTS[index]].StartingGess(x, i);
					auto st = chrono::high_resolution_clock::now();
					sequential::test[TESTS[index]].valgrad(x, g, i);
					auto diff = chrono::high_resolution_clock::now() - st;
					auto time = chrono::duration_cast<chrono::milliseconds>(diff);
					sequentialValgrads[r] = time.count();
					free(x);
					free(g);
				}
				//parallel valgrad
				{
					x = (double *)malloc(i * sizeof(double));
					g = (double *)malloc(i * sizeof(double));
					parallel::paralleltest[TESTS[index]].StartingGess(x, i);
					auto st = chrono::high_resolution_clock::now();
					parallel::paralleltest[TESTS[index]].valgrad(x, g, i);
					auto diff = chrono::high_resolution_clock::now() - st;
					auto time = chrono::duration_cast<chrono::milliseconds>(diff);
					parallelValgrads[r] = time.count();
					free(x);
					free(g);
				}
			}
			iSort(sequentialValues, REPETITIONS);
			iSort(sequentialGrads, REPETITIONS);
			iSort(sequentialValgrads, REPETITIONS);
			iSort(parallelValues, REPETITIONS);
			iSort(parallelGrads, REPETITIONS);
			iSort(parallelValgrads, REPETITIONS);

			cout << "n = " << i << endl;
			cout << "Sequential Value: " << sequentialValues[REPETITIONS / 2] << endl;
			cout << "Parallel Value: " << parallelValues[REPETITIONS / 2] << endl;
			cout << "Sequential Grad: " << sequentialGrads[REPETITIONS / 2] << endl;
			cout << "Parallel Grad: " << parallelGrads[REPETITIONS / 2] << endl;
			cout << "Sequential Valgrad: " << sequentialValgrads[REPETITIONS / 2] << endl;
			cout << "Parallel Valgrad: " << parallelValgrads[REPETITIONS / 2] << endl;

			cout << endl;

			free(sequentialValues);
			free(sequentialGrads);
			free(sequentialValgrads);
			free(parallelValues);
			free(parallelGrads);
			free(parallelValgrads);
		}
	}


}

int main(void)
{	
		runPerformanceTests();
}