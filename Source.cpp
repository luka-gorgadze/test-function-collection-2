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

vector<double> results;

#include "genericFunctions.h"
#include "testFunctions.h"
#include "parallelTestFunctions.h"

void run_cg_descent_Tests(struct testFunction* t, INT number)
{

	int repNumber = 10;
	double running_times[10];
	cg_stats Stats;
	char command;

	puts("Sequential and parallel implementations of test functions are available. To run parallel version enter P, To run sequential version enter S");
	scanf("%c", &command);
	if (command == 'S') {
		makeTestCollection(test);
	}
	if (command == 'P') {
		makeParallelTestCollection(test);
	}
	if (command != 'S' && command != 'P') {
		return;
	}

	puts("The following tests are available. plz enter test number");
	INT i;
	for (i = 0; i < number; i++)
	{
		if (test[i].size != 0) {
			printf("%d. %s, size of test: %d\n", i + 1, test[i].name, test[i].size);
		}
	}

	scanf("%d", &i);
	getchar();
	i--;
	n = test[i].size;
	printf("Chosen: %s, size of test: %d\n ", test[i].name, test[i].size);

	double *x;
	x = (double *)malloc(n * sizeof(double));
	int j;
	for (j = 0; j < repNumber; j++)
	{
		test[i].StartingGess(x, n);
		cg_parameter Parm;
		cg_default(&Parm);

		//Run cg descent algorithm and save its running time
		auto st = chrono::high_resolution_clock::now();
		cg_descent(x, n, &Stats, &Parm, 1.e-6, test[i].value, test[i].grad, test[i].valgrad, NULL);
		auto diff = chrono::high_resolution_clock::now() - st;
		auto time = chrono::duration_cast<chrono::milliseconds>(diff);
		running_times[j] = time.count();
	}
	iSort(running_times, 10);

	//Print results to file
	freopen("results.txt", "w", stdout);
	{
		printf("Median time:   %f\n", running_times[repNumber / 2]);
		printf("Functions value:   %f\n", Stats.f);
		int quarter = n / 4;
		int ii;
		for (ii = 0; ii < quarter; ii++)
		{
			printf("x[%3d] = %-15.11f     x[%3d] = %15.11f     x[%3d] = %15.11f     x[%3d] = %2.11f  \n",
				ii, x[ii], quarter + ii, x[quarter + ii], 2 * quarter + ii, x[2 * quarter + ii], 3 * quarter + ii, x[3 * quarter + ii]);
		}
	}
	free(x);
}


int main(void)
{
	run_cg_descent_Tests(test, TestNumber);
}