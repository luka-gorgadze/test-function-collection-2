#pragma once
INT n;

void iSort(double* a, int size)
{
	if (size > 1)
	{
		int i, j;
		double tmp;
		for (i = 1; i < size; ++i)
		{
			j = i - 1;
			tmp = a[i];
			while (j > -1 && a[j] > tmp)
			{
				a[j + 1] = a[j];
				--j;
			}
			a[j + 1] = tmp;
		}
	}
}

struct testFunction
{
	int size;
	char name[40];
	double(*value)
		(
			double   *x,
			INT       n
			);
	void(*grad)
		(
			double    *g,
			double    *x,
			INT        n
			);
	double(*valgrad)
		(
			double    *g,
			double    *x,
			INT        n
			);
	void(*StartingGess)
		(
			double   *x,
			INT       n
			);
};