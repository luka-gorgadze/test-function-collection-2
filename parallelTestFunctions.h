#pragma once

namespace parallel {

	//========== ARWHEAD ===========
	double ARWHEADvalgradRanged(double *g, double* x, double* sum, INT index, INT begin, INT end)
	{
		INT i;
		double fx = 0.0;
		double group1 = 0.0;
		for (i = begin; i < end; i++)
		{
			group1 = pow(x[i], 2) + pow(x[index], 2);
			fx += group1*group1 + 3.0 - 4.0*x[i];
			g[i] = -4.0 + 4.0*group1*x[i];
			*sum += 4.0*group1*x[index];
		}
		return fx;
	}
	double ARWHEADvalgrad(double *g, double *x, INT n)
	{
		g[n - 1] = 0;
		return valgradParallelOneCongestion(ARWHEADvalgradRanged, g, x, n - 1, n - 1);
	}

	void ARWHEADgradRanged(double *g, double *x, double *sum, INT index, INT begin, INT end)
	{
		INT i;
		double group1 = 0.0;

		for (i = begin; i < end; i++)
		{
			group1 = pow(x[i], 2) + pow(x[index], 2);
			g[i] = -4.0 + 4.0*group1*x[i];
			*sum += 4.0*group1*x[index];
		}
	}
	void ARWHEADgrad(double *g, double *x, INT n)
	{
		g[n - 1] = 0;
		gradParallelOneCongestion(ARWHEADgradRanged, g, x, n - 1, n - 1);
	}

	double ARWHEADvalueRanged(double *x, INT index, INT begin, INT end)
	{
		INT i;
		double fx = 0.0;
		double group1 = 0.0;
		for (i = begin; i < end; ++i)
		{
			group1 = pow(x[i], 2) + pow(x[index], 2);
			fx += group1*group1 + 3.0 - 4.0*x[i];
		}
		return fx;
	}
	double ARWHEADvalue(double *x, INT n)
	{
		return valueParallelOneCongestion(ARWHEADvalueRanged, x, n - 1, n - 1);
	}

	void ARWHEADStartingGess(double *x, INT n)
	{
		initializer(n - 1, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 1.0;
	}

	//========== BDQRTIC ===========
	double BDQRTICvalgradRanged(double *g, double* x, double* sum, INT index, INT begin, INT end)
	{
		INT i;
		double fx = 0.0;
		double group1 = 0.0, group2 = 0.0;
		for (i = begin; i < end; i++)
			g[i] = 0;
		for (i = begin; i < end - 4; i++)
		{
			group1 = 3.0 - 4.0*x[i];
			group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
				+ 4 * pow(x[i + 3], 2) + 5 * pow(x[index], 2);
			fx += group1*group1 + group2*group2;
			g[i] += -8.0*group1 + 4.0*group2*x[i];
			g[i + 1] += 8.0*group2*x[i + 1];
			g[i + 2] += 12.0*group2*x[i + 2];
			g[i + 3] += 16.0*group2*x[i + 3];
			*sum += 20.0*group2*x[index];
		}
		return fx;
	}
	double BDQRTICvalgrad(double *g, double *x, INT n)
	{
		double fx = valgradParallelOneCongestion(BDQRTICvalgradRanged, g, x, n, n - 1);
		INT i = block_size - 4;
		double group1 = 0.0, group2 = 0.0;
		for (int k = 1; k < num_threads; ++k) {
			for (int j = 0; j < 4; ++j) {
				group1 = 3.0 - 4.0*x[i];
				group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
					+ 4 * pow(x[i + 3], 2) + 5 * pow(x[n - 1], 2);
				g[i] += -8.0*group1 + 4.0*group2*x[i];
				g[i + 1] += 8.0*group2*x[i + 1];
				g[i + 2] += 12.0*group2*x[i + 2];
				g[i + 3] += 16.0*group2*x[i + 3];
				++i;
			}
			i += (block_size - 4);
		}
		return fx;
	}

	void BDQRTICgradRanged(double *g, double *x, double *sum, INT index, INT begin, INT end)
	{

		INT i;
		double group1 = 0.0, group2 = 0.0;
		for (i = begin; i < end; i++)
			g[i] = 0;
		for (i = begin; i < end - 4; i++)
		{
			group1 = 3.0 - 4.0*x[i];
			group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
				+ 4 * pow(x[i + 3], 2) + 5 * pow(x[index], 2);
			g[i] += -8.0*group1 + 4.0*group2*x[i];
			g[i + 1] += 8.0*group2*x[i + 1];
			g[i + 2] += 12.0*group2*x[i + 2];
			g[i + 3] += 16.0*group2*x[i + 3];
			*sum += 20.0*group2*x[index];
		}
	}
	void BDQRTICgrad(double *g, double *x, INT n)
	{
		gradParallelOneCongestion(BDQRTICgradRanged, g, x, n, n - 1);
		INT i = block_size - 4;
		double group1 = 0.0, group2 = 0.0;
		for (int k = 1; k < num_threads; ++k) {
			for (int j = 0; j < 4; ++j) {
				group1 = 3.0 - 4.0*x[i];
				group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
					+ 4 * pow(x[i + 3], 2) + 5 * pow(x[n - 1], 2);
				g[i] += -8.0*group1 + 4.0*group2*x[i];
				g[i + 1] += 8.0*group2*x[i + 1];
				g[i + 2] += 12.0*group2*x[i + 2];
				g[i + 3] += 16.0*group2*x[i + 3];
				++i;
			}
			i += (block_size - 4);
		}
	}

	double BDQRTICvalueRanged(double *x, INT index, INT begin, INT end)
	{
		INT i;
		double fx = 0.0;
		double group1 = 0.0, group2 = 0.0;
		for (i = begin; i < end - 4; i++)
		{
			group1 = 3.0 - 4.0*x[i];
			group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
				+ 4 * pow(x[i + 3], 2) + 5 * pow(x[index - 1], 2);
			fx += group1*group1 + group2*group2;
		}
		return fx;
	}
	double BDQRTICvalue(double *x, INT n)
	{
		return valueParallelOneCongestion(BDQRTICvalueRanged, x, n - 1, n);
	}

	void BDQRTICStartingGess(double *x, INT n)
	{
		initializer(n - 1, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 1.0;
	}

	//========== COSINE ===========
	double COSINEvalgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		double fx = 0.0;

		for (i = begin; i < end; i++)
			g[i] = 0.0;

		for (i = begin; i < end - 1; i++)
		{
			item = -0.5*x[i + 1] + x[i] * x[i];
			fx += cos(item);
			g[i] -= 2.0*sin(item)*x[i];
			g[i + 1] += 0.5*sin(item);
		}
		return fx;
	}
	double COSINEvalgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		double fx = valgradParallel(COSINEvalgradRanged, g, x, n);
		double item;
		for (int i = block_size - 1; i < n - 1; i = i + block_size) {
			item = -0.5*x[i + 1] + x[i] * x[i];
			g[i] -= 2.0*sin(item)*x[i];
			g[i + 1] += 0.5*sin(item);
			fx += cos(item);
		}
		return fx;
	}

	void COSINEgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		for (i = begin; i < end; i++)
			g[i] = 0.0;

		for (i = begin; i < end - 1; i++)
		{
			item = -0.5*x[i + 1] + x[i] * x[i];
			g[i] -= 2.0*sin(item)*x[i];
			g[i + 1] += 0.5*sin(item);
		}
	}
	void COSINEgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		gradParallel(COSINEgradRanged, g, x, n);
		double item;
		for (int i = block_size - 1; i < n - 1; i = i + block_size) {
			item = -0.5*x[i + 1] + x[i] * x[i];
			g[i] -= 2.0*sin(item)*x[i];
			g[i + 1] += 0.5*sin(item);
		}
	}

	double COSINEvalueRanged
	(
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		double fx = 0.0;

		for (i = begin; i < end; i++)
		{
			item = -0.5*x[i + 1] + x[i] * x[i];
			fx += cos(item);
		}
		return fx;
	}
	double COSINEvalue
	(
		double *x,
		INT n
	)
	{
		return valueParallel(COSINEvalueRanged, x, n - 1);
	}

	void COSINEStartingGess
	(
		double *x,
		INT n
	)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 1.0;
	}

	//========== DIXON3DQ ===========
	double DIXON3DQvalgradRanged(double *g, double* x, INT begin, INT end)
	{
		INT i;
		double item;
		double fx = (0.0);
		for (i = begin + 1; i <= end; i++)
			g[i] = 0.0;
		for (i = begin + 1; i < end; i++)
		{
			item = (x[i] - x[i + 1]);
			g[i] += 2.0*item;
			g[i + 1] -= 2.0*item;
			item = (x[i - 1] - x[i]);
			fx += item*item;
		}
		return fx;

	}
	double DIXON3DQvalgrad(double *g, double *x, INT n)
	{
		if (num_threads == 1)
		{
			INT i;
			double item;
			double fx = (0.0);
			item = x[0] - 1;
			fx += item*item;
			g[0] = 2.0*item;
			for (i = 1; i < n - 1; i++)
				g[i] = 0.0;
			item = x[n - 1] - 1;
			fx += item*item;
			g[n - 1] = 2.0*item;
			for (i = 1; i < n - 1; i++)
			{
				item = (x[i] - x[i + 1]);
				fx += item*item;
				g[i] += 2.0*item;
				g[i + 1] -= 2.0*item;
			}
			return fx;
		}
		double fx = 0.0;
		double item;
		item = x[0] - 1;
		fx += item*item;
		g[0] = 2.0*item;
		item = x[n - 1] - 1;
		fx += item*item;
		fx += valgradParallel(DIXON3DQvalgradRanged, g, x, n - 1);
		INT block_start = 0, i;
		for (i = 0; i < (num_threads - 1); ++i)
		{
			block_start += block_size;
			item = (x[block_start] - x[block_start + 1]);
			g[block_start] += 2.0*item;
			g[block_start + 1] -= 2.0*item;
		}
		item = x[n - 1] - 1;
		g[n - 1] += 2.0*item;
		return fx;
	}

	void DIXON3DQgradRanged(double *g, double *x, INT begin, INT end)
	{
		INT i;
		double item;
		for (i = begin + 1; i <= end; i++)
			g[i] = 0.0;
		for (i = begin + 1; i < end; i++)
		{
			item = (x[i] - x[i + 1]);
			g[i] += 2.0*item;
			g[i + 1] -= 2.0*item;
		}
	}
	void DIXON3DQgrad(double *g, double *x, INT n)
	{
		double item;
		item = x[0] - 1;
		g[0] = 2.0*item;
		gradParallel(DIXON3DQgradRanged, g, x, n - 1);
		INT block_start = 0, i;
		for (i = 0; i < (num_threads - 1); ++i)
		{
			block_start += block_size;
			item = (x[block_start] - x[block_start + 1]);
			g[block_start] += 2.0*item;
			g[block_start + 1] -= 2.0*item;
		}

		item = x[n - 1] - 1;
		g[n - 1] += 2.0*item;
	}

	double DIXON3DQvalueRanged(double *x, INT begin, INT end)
	{
		INT i;
		double item;
		double fx = (0.0);
		for (i = begin; i < end; i++)
		{
			item = (x[i] - x[i + 1]);
			fx += item*item;
		}
		return fx;
	}
	double DIXON3DQvalue(double *x, INT n)
	{
		if (num_threads == 1)
		{
			INT i;
			double item;
			double fx = (0.0);
			item = x[0] - 1;
			fx += item*item;
			item = x[n - 1] - 1;
			fx += item*item;
			for (i = 1; i < n - 1; i++)
			{
				item = (x[i] - x[i + 1]);
				fx += item*item;
			}
			return fx;
		}
		double fx = 0.0;
		double item;
		item = x[0] - 1;
		fx += item*item;
		item = x[n - 1] - 1;
		fx += item*item;
		fx += valueParallel(DIXON3DQvalueRanged, x, n - 1);
		return fx;
	}

	void DIXON3DQStartingGess(double *x, INT n)
	{
		initializer(n - 1, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = -1.0;
	}

	//========== DQDRTIC ===========
	double DQDRTICvalgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double fx = 0.0;
		for (i = begin; i < end; i++)
			g[i] = 0;
		for (i = begin; i < end - 2; i++)
		{
			fx += pow(x[i], 2) + 100 * pow(x[i + 1], 2) + 100 * pow(x[i + 2], 2);
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
		}
		return fx;
	}
	double DQDRTICvalgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		double fx = valgradParallel(DQDRTICvalgradRanged, g, x, n);
		for (int i = block_size - 2; i < n - 2; i = i + block_size) {
			fx += pow(x[i], 2) + 100 * pow(x[i + 1], 2) + 100 * pow(x[i + 2], 2);
			fx += pow(x[i + 1], 2) + 100 * pow(x[i + 2], 2) + 100 * pow(x[i + 3], 2);
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
			g[i + 1] += 2 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
			g[i + 3] += 200 * x[i + 3];
		}
		return fx;
	}

	void DQDRTICgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		for (i = begin; i < end; i++)
			g[i] = 0;
		for (i = begin; i < end - 2; i++)
		{
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
		}
	}
	void DQDRTICgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		gradParallel(DQDRTICgradRanged, g, x, n);
		for (int i = block_size - 2; i < n - 2; i = i + block_size) {
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
			g[i + 1] += 2 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
			g[i + 3] += 200 * x[i + 3];
		}
	}

	double DQDRTICvalueRanged
	(
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double fx = 0.0;
		for (i = begin; i < end; i++)
		{
			fx += pow(x[i], 2) + 100 * pow(x[i + 1], 2) + 100 * pow(x[i + 2], 2);
		}
		return fx;
	}
	double DQDRTICvalue
	(
		double *x,
		INT n
	)
	{
		return valueParallel(DQDRTICvalueRanged, x, n - 2);
	}

	void DQDRTICStartingGess
	(
		double *x,
		INT n
	)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 3.0;
	}

	//========== DQRTIC ===========
	double DQRTICvalgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double f = 0.0;
		double item;
		for (i = begin; i < end; i++)
		{
			item = x[i] - i - 1;
			f += pow(item, 4);
			g[i] = 4 * pow(item, 3);
		}
		return f;
	}
	double DQRTICvalgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		return valgradParallel(DQRTICvalgradRanged, g, x, n);
	}

	void DQRTICgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		for (i = begin; i < end; i++)
		{
			item = x[i] - i - 1;
			g[i] = 4 * pow(item, 3);
		}
	}
	void DQRTICgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		gradParallel(DQRTICgradRanged, g, x, n);
	}

	double DQRTICvalueRanged
	(
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double f = 0.0;
		double item;
		for (i = begin; i < end; i++)
		{
			item = x[i] - i - 1;
			f += pow(item, 4);
		}
		return f;
	}
	double DQRTICvalue
	(
		double *x,
		INT n
	)
	{
		return valueParallel(DQRTICvalueRanged, x, n);
	}

	void DQRTICStartingGess
	(
		double *x,
		INT n
	)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 2.0;
	}

	//========== EDENSCH ===========
	double EDENSCHvalgradRanged(double*g, double*x, INT begin, INT end)
	{
		double item1, item2, item3;
		for (int i = begin; i < end; ++i)
			g[i] = 0;
		double fx = 0.0;
		for (int i = begin; i < end - 1; ++i)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			fx += 16 + pow(item1, 4) + pow(item2, 2) + pow(item3, 2);
			g[i] += 4 * pow(item1, 3) + 2 * item2*x[i + 1];
			g[i + 1] += 2 * item2*(x[i] - 2.0) + 2 * item3;
		}
		return fx;
	}
	double EDENSCHvalgrad(double *g, double*x, INT n)
	{
		double fx = valgradParallel(EDENSCHvalgradRanged, g, x, n);
		double item1, item2, item3;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			g[i] += 4 * pow(item1, 3) + 2 * item2*x[i + 1];
			g[i + 1] += 2 * item2*(x[i] - 2.0) + 2 * item3;
			fx += 16 + pow(item1, 4) + pow(item2, 2) + pow(item3, 2);
		}
		return fx;
	}

	void EDENSCHgradRanged(double*g, double*x, INT begin, INT end)
	{
		INT i;
		for (i = begin; i < end; i++)
			g[i] = 0;
		double item1, item2, item3;
		for (i = begin; i < end - 1; i++)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			g[i] += 4 * pow(item1, 3) + 2 * item2*x[i + 1];
			g[i + 1] += 2 * item2*(x[i] - 2.0) + 2 * item3;
		}
	}
	void EDENSCHgrad(double*g, double*x, INT n)
	{
		gradParallel(EDENSCHgradRanged, g, x, n);
		double item1, item2, item3;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			g[i] += 4 * pow(item1, 3) + 2 * item2*x[i + 1];
			g[i + 1] += 2 * item2*(x[i] - 2.0) + 2 * item3;
		}
	}

	double EDENSCHvalueRanged(double *x, INT begin, INT end)
	{
		INT i;
		double fx = 0.0;
		double item1, item2, item3;
		for (i = begin; i < end; i++)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			fx += 16 + pow(item1, 4) + pow(item2, 2) + pow(item3, 2);
		}
		return fx;
	}
	double EDENSCHvalue(double *x, INT n)
	{
		return 	valueParallel(EDENSCHvalueRanged, x, n - 1);

	}

	void EDENSCHStartingGess(double* x, INT n)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 0.0;
	}

	//========== EG2 ===========
	double EG2valgradRanged
	(
		double *g,
		double *x,
		double *sum,
		INT index,
		INT begin,
		INT end
	)
	{
		INT i;
		for (i = begin; i < end; i++)
			g[i] = 0;
		double fx = 0;
		double item;
		for (i = begin; i < end; i++)
		{
			item = x[index] + pow(x[i], 2) - 1;;
			fx += sin(item);
			*sum += cos(item);
			g[i] += 2 * cos(item)*x[i];
		}
		return fx;
	}
	double EG2valgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		g[n - 1] = cos(pow(x[n - 1], 2))*x[n - 1];
		return 0.5*sin(pow(x[n - 1], 2)) + valgradParallelOneCongestion(EG2valgradRanged, g, x, n - 1, 0);
	}

	void EG2gradRanged
	(
		double *g,
		double *x,
		double *sum,
		INT index,
		INT begin,
		INT end
	)
	{
		INT i;
		for (i = begin; i < end; i++)
			g[i] = 0;
		double item;
		for (i = begin; i < end; i++)
		{
			item = x[index] + pow(x[i], 2) - 1;;
			*sum += cos(item);
			g[i] += 2 * cos(item)*x[i];
		}
	}
	void EG2grad
	(
		double *g,
		double *x,
		INT n
	)
	{
		g[n - 1] = cos(pow(x[n - 1], 2))*x[n - 1];
		gradParallelOneCongestion(EG2gradRanged, g, x, n - 1, 0);
	}

	double EG2valueRanged
	(
		double *x,
		INT index,
		INT begin,
		INT end
	)
	{
		INT i;
		double fx = 0;
		double item;
		for (i = begin; i < end; i++)
		{
			item = x[index] + pow(x[i], 2) - 1;;
			fx += sin(item);
		}
		return fx;
	}
	double EG2value
	(
		double *x,
		INT n
	)
	{
		return 0.5*sin(pow(x[n - 1], 2)) + valueParallelOneCongestion(EG2valueRanged, x, n - 1, 0);
	}

	void EG2StartingGess
	(
		double *x,
		INT n
	)
	{
		initializer(2000, n - 1);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 0.;
	}

	//========== ENGVAL1 ===========
	double ENGVAL1valgradRanged(double*g, double*x, INT begin, INT end)
	{
		double item;
		for (int i = begin; i < end; ++i)
			g[i] = 0;
		double fx = 0.0;
		for (int i = begin; i < end - 1; ++i)
		{
			item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
			fx += item*item + (3 - 4.0*x[i]);
			g[i] += 4.0*item*x[i] - 4.0;
			g[i + 1] += 4.0*item*x[i + 1];
		}
		return fx;
	}
	double ENGVAL1valgrad(double *g, double*x, INT n)
	{
		double fx = valgradParallel(ENGVAL1valgradRanged, g, x, n);
		double item;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
			fx += item*item + (3 - 4.0*x[i]);
			g[i] += 4.0*item*x[i] - 4.0;
			g[i + 1] += 4.0*item*x[i + 1];
		}
		return fx;
	}

	void ENGVAL1gradRanged(double*g, double*x, INT begin, INT end)
	{
		INT i;
		for (i = begin; i < end; i++)
			g[i] = 0.0;
		double item;
		for (i = begin; i < end - 1; i++)
		{
			item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
			g[i] += 4.0*item*x[i] - 4.0;
			g[i + 1] += 4.0*item*x[i + 1];
		}
	}
	void ENGVAL1grad(double*g, double*x, INT n)
	{
		gradParallel(ENGVAL1gradRanged, g, x, n);
		double item;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
			g[i] += 4.0*item*x[i] - 4.0;
			g[i + 1] += 4.0*item*x[i + 1];
		}
	}

	double ENGVAL1valueRanged(double *x, INT begin, INT end)
	{
		INT i;
		double item;
		double fx = (0.0);
		for (i = begin; i < end; i++)
		{
			item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
			fx += item*item + (3 - 4.0*x[i]);
		}
		return fx;

	}
	double ENGVAL1value(double *x, INT n)
	{
		return 	valueParallel(ENGVAL1valueRanged, x, n - 1);
	}

	void ENGVAL1StartingGess(double *x, INT n)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 2.0;
	}

	//========== FLETCHER ===========
	double FLETCHRvalgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		double fx = 0.0;

		for (i = begin; i < end; i++)
			g[i] = 0.0;

		for (i = begin; i < end - 1; i++)
		{
			item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
			fx += item*item;
			g[i] += 20.0*item*(-2.0*x[i] - 1.0);
			g[i + 1] += 20.0*item;
		}
		return fx;
	}
	double FLETCHRvalgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		double fx = valgradParallel(FLETCHRvalgradRanged, g, x, n);
		double item;
		for (int i = block_size - 1; i < n - 1; i = i + block_size) {
			item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
			fx += item*item;
			g[i] += 20.0*item*(-2.0*x[i] - 1.0);
			g[i + 1] += 20.0*item;
		}
		return 100.0*fx;
	}

	void FLETCHRgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		for (i = begin; i < end; i++)
			g[i] = 0.0;

		for (i = begin; i < end - 1; i++)
		{
			item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
			g[i] += 20.0*item*(-2.0*x[i] - 1.0);
			g[i + 1] += 20.0*item;
		}
	}
	void FLETCHRgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		gradParallel(FLETCHRgradRanged, g, x, n);
		double item;
		for (int i = block_size - 1; i < n - 1; i = i + block_size)
		{
			item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
			g[i] += 20.0*item*(-2.0*x[i] - 1.0);
			g[i + 1] += 20.0*item;
		}
	}

	double FLETCHRvalueRanged
	(
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		double fx = 0.0;

		for (i = begin; i < end; i++)
		{
			item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
			fx += item*item;
		}
		return 100.0*fx;
	}
	double FLETCHRvalue
	(
		double *x,
		INT n
	)
	{
		return valueParallel(FLETCHRvalueRanged, x, n - 1);
	}

	void FLETCHRStartingGess
	(
		double *x,
		INT n
	)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 1.0;
	}

	//========== FREUROTH ===========
	double FREUROTHvalgradRanged(double*g, double*x, INT begin, INT end)
	{
		double item1, item2;
		for (int i = begin; i < end; ++i)
			g[i] = 0;
		double fx = 0.0;
		for (int i = begin; i < end - 1; ++i)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
			fx += item1*item1 + item2*item2;
			g[i] += 2.0*item1 + 2.0*item2;
			g[i + 1] += 2.0*item1*(10 * x[i + 1] - 3.0*x[i + 1] * x[i + 1] - 2.0) +
				2.0*item2*(2 * x[i + 1] + 3.0*x[i + 1] * x[i + 1] - 14.0);
		}
		return fx;
	}
	double FREUROTHvalgrad(double *g, double*x, INT n)
	{
		double fx = valgradParallel(FREUROTHvalgradRanged, g, x, n);
		double item1, item2;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
			fx += item1*item1 + item2*item2;
			g[i] += 2.0*item1 + 2.0*item2;
			g[i + 1] += 2.0*item1*(10 * x[i + 1] - 3.0*x[i + 1] * x[i + 1] - 2.0) +
				2.0*item2*(2 * x[i + 1] + 3.0*x[i + 1] * x[i + 1] - 14.0);
		}
		return fx;
	}

	void FREUROTHgradRanged(double*g, double*x, INT begin, INT end)
	{
		INT i;
		for (i = begin; i < end; i++)
			g[i] = 0.0;
		double item1, item2;
		for (i = begin; i < end - 1; i++)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
			g[i] += 2.0*item1 + 2.0*item2;
			g[i + 1] += 2.0*item1*(10 * x[i + 1] - 3.0*x[i + 1] * x[i + 1] - 2.0) +
				2.0*item2*(2 * x[i + 1] + 3.0*x[i + 1] * x[i + 1] - 14.0);
		}
	}
	void FREUROTHgrad(double*g, double*x, INT n)
	{
		gradParallel(FREUROTHgradRanged, g, x, n);
		double item1, item2;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
			g[i] += 2.0*item1 + 2.0*item2;
			g[i + 1] += 2.0*item1*(10 * x[i + 1] - 3.0*x[i + 1] * x[i + 1] - 2.0) +
				2.0*item2*(2 * x[i + 1] + 3.0*x[i + 1] * x[i + 1] - 14.0);
		}
	}

	double FREUROTHvalueRanged(double *x, INT begin, INT end)
	{
		INT i;
		double item1, item2;
		double fx = (0.0);
		for (i = begin; i < end; i++)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
			fx += item1*item1 + item2*item2;
		}
		return fx;

	}
	double FREUROTHvalue(double *x, INT n)
	{
		return 	valueParallel(FREUROTHvalueRanged, x, n - 1);

	}

	void FREUROTHStartingGess(double *x, INT n)
	{
		initializer(n, 2000);
		INT i;
		x[0] = 0.5;
		x[1] = -2.0;
		for (i = 2; i < n; i++)
			x[i] = 0.0;
	}

	//========== GENHUMPS ===========
	double GENHUMPSvalgradRanged(double*g, double*x, INT begin, INT end)
	{
		double item1, item2;
		for (int i = begin; i < end; ++i)
			g[i] = 0;
		double fx = 0.0;
		for (int i = begin; i < end - 1; ++i)
		{
			item1 = sin(2.0*x[i]);
			item2 = sin(2.0*x[i + 1]);
			fx += item1*item1 *item2*item2 + 0.05*(pow(x[i], 2) + pow(x[i + 1], 2));
			g[i] += 4.0*item1*cos(2.0*x[i]) *item2*item2 + 0.1*x[i];
			g[i + 1] += 4.0*item2*cos(2.0*x[i + 1])* item1*item1 + 0.1*x[i + 1];
		}
		return fx;
	}

	double GENHUMPSvalgrad(double *g, double*x, INT n)
	{
		double fx = valgradParallel(GENHUMPSvalgradRanged, g, x, n);
		double item1, item2;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = sin(2.0*x[i]);
			item2 = sin(2.0*x[i + 1]);
			fx += item1*item1 *item2*item2 + 0.05*(pow(x[i], 2) + pow(x[i + 1], 2));
			g[i] += 4.0*item1*cos(2.0*x[i]) *item2*item2 + 0.1*x[i];
			g[i + 1] += 4.0*item2*cos(2.0*x[i + 1])* item1*item1 + 0.1*x[i + 1];
		}
		return fx;
	}

	void GENHUMPSgradRanged(double*g, double*x, INT begin, INT end)
	{
		INT i;
		for (i = begin; i < end; i++)
			g[i] = 0.0;
		double item1, item2;
		for (i = begin; i < end - 1; i++)
		{
			item1 = sin(2.0*x[i]);
			item2 = sin(2.0*x[i + 1]);
			g[i] += 4.0*item1*cos(2.0*x[i]) *item2*item2 + 0.1*x[i];
			g[i + 1] += 4.0*item2*cos(2.0*x[i + 1])* item1*item1 + 0.1*x[i + 1];
		}
	}
	void GENHUMPSgrad(double*g, double*x, INT n)
	{
		gradParallel(GENHUMPSgradRanged, g, x, n);
		double item1, item2;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = sin(2.0*x[i]);
			item2 = sin(2.0*x[i + 1]);
			g[i] += 4.0*item1*cos(2.0*x[i]) *item2*item2 + 0.1*x[i];
			g[i + 1] += 4.0*item2*cos(2.0*x[i + 1])* item1*item1 + 0.1*x[i + 1];
		}
	}

	double GENHUMPSvalueRanged(double *x, INT begin, INT end)
	{
		INT i;
		double item1, item2;
		double fx = (0.0);
		for (i = begin; i < end; i++)
		{
			item1 = sin(2.0*x[i]);
			item2 = sin(2.0*x[i + 1]);
			fx += item1*item1 *item2*item2 + 0.05*(pow(x[i], 2) + pow(x[i + 1], 2));
		}
		return fx;

	}
	double GENHUMPSvalue(double *x, INT n)
	{
		return 	valueParallel(GENHUMPSvalueRanged, x, n - 1);

	}

	void GENHUMPSStartingGess(double *x, INT n)
	{
		initializer(n, 2000);
		INT i;
		x[0] = -506.0;
		for (i = 1; i < n; i++)
			x[i] = 506.2;
	}

	//========== GENROSE ===========
	double GENROSEvalgradRanged(double *g, double* x, INT begin, INT end)
	{
		INT i;
		double fx = 1.0;
		double item1, item2;
		for (i = begin; i < end; i++)
			g[i] = 0;
		for (i = begin + 1; i < end; i++)
		{
			item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
			item2 = x[i] - 1.0;
			fx += item1*item1 + item2*item2;
			g[i - 1] -= 40.0*item1 * x[i - 1];
			g[i] += 20 * item1 + 2 * item2;
		}
		return fx;

	}
	double GENROSEvalgrad(double *g, double *x, INT n)
	{
		double fx = 0.0;
		double item1, item2;
		fx += valgradParallel(GENROSEvalgradRanged, g, x, n);
		INT block_start = 0, i;
		for (i = 0; i < (num_threads - 1); ++i)
		{
			block_start += block_size;
			item1 = 10.0 * (x[block_start] - x[block_start - 1] * x[block_start - 1]);
			item2 = x[block_start] - 1.0;
			g[block_start - 1] -= 40.0*item1 * x[block_start - 1];
			g[block_start] += 20 * item1 + 2 * item2;
		}
		return fx;
	}

	void GENROSEgradRanged(double *g, double *x, INT begin, INT end)
	{
		INT i;
		double item1, item2;
		for (i = begin; i < end; i++)
			g[i] = 0;
		for (i = begin + 1; i < end; i++)
		{
			item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
			item2 = x[i] - 1.0;
			g[i - 1] -= 40.0*item1 * x[i - 1];
			g[i] += 20 * item1 + 2 * item2;
		}

	}
	void GENROSEgrad(double *g, double *x, INT n)
	{
		double item1, item2;
		gradParallel(GENROSEgradRanged, g, x, n);
		INT block_start = 0, i;
		for (i = 0; i < (num_threads - 1); ++i)
		{
			block_start += block_size;
			item1 = 10.0 * (x[block_start] - x[block_start - 1] * x[block_start - 1]);
			item2 = x[block_start] - 1.0;
			g[block_start - 1] -= 40.0*item1 * x[block_start - 1];
			g[block_start] += 20 * item1 + 2 * item2;
		}
	}

	double GENROSEvalueRanged(double *x, INT begin, INT end)
	{
		INT i;
		double fx = 1.0;
		double item1, item2;
		for (i = begin + 1; i < end; i++)
		{
			item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
			item2 = x[i] - 1.0;
			fx += item1*item1 + item2*item2;
		}
		return fx;
	}
	double GENROSEvalue(double *x, INT n)
	{
		return valueParallel(GENROSEvalueRanged, x, n);
	}

	void GENROSEStartingGess(double *x, INT n)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 1.0 / n + 1;

	}

	//========== LIARWHD ===========
	double LIARWDHvalgradRanged(double *g, double* x, double* sum, INT index, INT begin, INT end)
	{
		INT i;
		double item1, item2;
		double fx = (0.0);
		for (i = begin; i < end; i++)
		{
			item1 = 2.0*(x[i] * x[i] - x[0]);
			item2 = (x[i] - 1);
			fx += item1*item1 + item2*item2;
			g[i] = 8.0*item1*x[i] + 2.0*item2;
			*sum -= 4.0*item1;
		}
		return fx;
	}
	double LIARWDHvalgrad(double *g, double *x, INT n)
	{
		return valgradParallelOneCongestion(LIARWDHvalgradRanged, g, x, n, 0);
	}

	void LIARWDHgradRanged(double *g, double *x, double *sum, INT index, INT begin, INT end)
	{
		INT i;
		double item1, item2;
		for (i = begin; i < end; i++)
		{
			item1 = 2.0*(x[i] * x[i] - x[index]);
			item2 = (x[i] - 1);
			g[i] = 8.0*item1*x[i] + 2.0*item2;
			*sum -= 4.0*item1;
		}
	}
	void LIARWDHgrad(double *g, double *x, INT n)
	{
		gradParallelOneCongestion(LIARWDHgradRanged, g, x, n, 0);
	}

	double LIARWDHvalueRanged(double *x, INT index, INT begin, INT end)
	{
		INT i;
		double item1, item2;
		double fx = (0.0);
		for (i = begin; i < end; i++)
		{
			item1 = 2.0*(x[i] * x[i] - x[0]);
			item2 = (x[i] - 1);
			fx += item1*item1 + item2*item2;
		}
		return fx;
	}
	double LIARWDHvalue(double *x, INT n)
	{
		return valueParallelOneCongestion(LIARWDHvalueRanged, x, n, 0);
	}

	void LIARWDHStartingGess(double *x, INT n)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i++)
			x[i] = 4.0;
	}

	//========== NONDQUAR ===========
	double NONDQUARvalgradRanged(double *g, double* x, double* sum, INT index, INT begin, INT end)
	{
		INT i;
		double item;
		double fx = (0.0);
		for (i = begin; i < end; i++)
			g[i] = 0.0;
		for (i = begin; i < end; i++)
		{
			item = x[i] + x[i + 1] + x[index];
			fx += pow(item, 4.0);
			g[i] += 4.0*pow(item, 3.0);
			g[i + 1] += 4.0*pow(item, 3.0);
			*sum += 4.0*pow(item, 3.0);
		}
		return fx;

	}
	double NONDQUARvalgrad(double *g, double *x, INT n)
	{
		double fx = 0.0;
		double item;
		g[n - 1] = 0;
		g[n - 2] = 0;
		fx += valgradParallelOneCongestion(NONDQUARvalgradRanged, g, x, n - 2, n - 1);
		item = x[0] - x[1];
		g[0] += 2.0*item;
		g[1] += 2.0*item;
		fx += item*item;
		item = x[n - 2] + x[n - 1];
		g[n - 2] += 2.0 * item;
		g[n - 1] += 2.0 * item;
		fx += item*item;
		return fx;

	}

	void NONDQUARgradRanged(double *g, double *x, double *sum, INT index, INT begin, INT end)
	{
		INT i;
		double item;
		for (i = begin; i < end; i++)
			g[i] = 0.0;
		for (i = begin; i < end; i++)
		{
			item = x[i] + x[i + 1] + x[index];
			g[i] += 4.0*pow(item, 3.0);
			g[i + 1] += 4.0*pow(item, 3.0);
			*sum += 4.0*pow(item, 3.0);

		}

	}
	void NONDQUARgrad(double *g, double *x, INT n)
	{
		int item;
		g[n - 1] = 0;
		g[n - 2] = 0;
		gradParallelOneCongestion(NONDQUARgradRanged, g, x, n - 2, n - 1);
		item = x[0] - x[1];
		g[0] += 2.0 * item;
		g[1] += -2.0 * item;
		item = x[n - 2] + x[n - 1];
		g[n - 2] += 2.0 * item;
		g[n - 1] += 2.0 * item;
	}

	double NONDQUARvalueRanged(double *x, INT index, INT begin, INT end)
	{
		INT i;
		double item;
		double fx = (0.0);
		for (i = begin; i < end; i++)
		{
			item = x[i] + x[i + 1] + x[index - 1];
			fx += pow(item, 4.0);
		}
		return fx;
	}
	double NONDQUARvalue(double *x, INT n)
	{
		if (num_threads == 1)
		{
			INT i;
			double item;
			double fx = (0.0);
			item = x[0] - x[1];
			fx += item*item;
			item = x[n - 2] + x[n - 1];
			fx += item*item;
			for (i = 0; i < n - 2; i++)
			{
				item = x[i] + x[i + 1] + x[n - 1];
				fx += pow(item, 4.0);
			}
			return fx;
		}
		double fx = 0.0;
		double item;
		fx += valueParallelOneCongestion(NONDQUARvalueRanged, x, n - 1, n);
		item = x[0] - x[1];
		fx += item*item;
		item = x[n - 2] + x[n - 1];
		fx += item*item;
		return fx;
	}

	void NONDQUARStartingGess(double *x, INT n)
	{
		initializer(n - 1, 2000);
		INT i;
		for (i = 0; i < n; i += 2)
		{
			x[i] = 1.0;
			x[i + 1] = -1.0;
		}
	}

	//========== POWER ===========
	double POWERvalgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		double fx = (0.0);

		for (i = begin; i < end; i++)
		{
			item = (i + 1)*x[i];
			fx += item*item;
			g[i] = 2.0*item*(i + 1);
		}
		return fx;
	}
	double POWERvalgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		return valgradParallel(POWERvalgradRanged, g, x, n);
	}

	void POWERgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		for (i = begin; i < end; i++)
		{
			item = (i + 1)*x[i];
			g[i] = 2.0*item*(i + 1);
		}
	}
	void POWERgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		gradParallel(POWERgradRanged, g, x, n);
	}

	double POWERvalueRanged
	(
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double item;
		double fx = (0.0);

		for (i = begin; i < end; i++)
		{
			item = (i + 1)*x[i];
			fx += item*item;
		}
		return fx;
	}
	double POWERvalue
	(
		double *x,
		INT n
	)
	{
		return valueParallel(POWERvalueRanged, x, n);
	}

	void POWERStartingGess
	(
		double *x,
		INT n
	)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i<n; i++)
			x[i] = 1.0;

	}

	//========== SROSENBR ===========
	double SROSENBRvalgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double fx = 0.0;
		for (i = begin; i < end; i += 2) {
			double t1 = 1.0 - x[i];
			double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
			g[i + 1] = 20.0 * t2;
			g[i] = -2.0 * (x[i] * g[i + 1] + t1);
			fx += t1 * t1 + t2 * t2;
		}
		return fx;
	}
	double SROSENBRvalgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		return valgradParallel(SROSENBRvalgradRanged, g, x, n);
	}

	void SROSENBRgradRanged
	(
		double *g,
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		for (i = begin; i < end; i += 2) {
			double t1 = 1.0 - x[i];
			double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
			g[i + 1] = 20.0 * t2;
			g[i] = -2.0 * (x[i] * g[i + 1] + t1);
		}
	}
	void SROSENBRgrad
	(
		double *g,
		double *x,
		INT n
	)
	{
		return gradParallel(SROSENBRgradRanged, g, x, n);
	}

	double SROSENBRvalueRanged
	(
		double *x,
		INT begin,
		INT end
	)
	{
		INT i;
		double fx = 0.0;
		for (i = begin; i < end; i += 2) {
			double t1 = 1.0 - x[i];
			double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
			fx += t1 * t1 + t2 * t2;
		}
		return fx;
	}
	double SROSENBRvalue
	(
		double *x,
		INT n
	)
	{
		return valueParallel(SROSENBRvalueRanged, x, n);
	}

	void SROSENBRStartingGess
	(
		double *x,
		INT n
	)
	{
		initializer(n, 2000);
		INT i;
		for (i = 0; i < n; i += 2)
		{
			x[i] = -1.2;
			x[i + 1] = 1.0;
		}
	}
}

const INT ParallelTestNumber = 45;
struct testFunction paralleltest[45];

void makeParallelTestCollection(struct testFunction* t)
{
	//1-ARWHEAD
	test[0].size = 5000;
	strcpy(test[0].name, "ARWHEAD");
	test[0].value = parallel::ARWHEADvalue;
	test[0].grad = parallel::ARWHEADgrad;
	test[0].valgrad = parallel::ARWHEADvalgrad;
	test[0].StartingGess = parallel::ARWHEADStartingGess;
	//2-BDQRTIC
	test[1].size = 5000;
	strcpy(test[1].name, "BDQRTIC");
	test[1].value = parallel::BDQRTICvalue;
	test[1].grad = parallel::BDQRTICgrad;
	test[1].valgrad = parallel::BDQRTICvalgrad;
	test[1].StartingGess = parallel::BDQRTICStartingGess;
	//3-BROYDN7D
	//4-BRYBND
	//5-CHAINWOO
	//6-COSINE
	test[5].size = 10000;
	strcpy(test[5].name, "COSINE");
	test[5].value = parallel::COSINEvalue;
	test[5].grad = parallel::COSINEgrad;
	test[5].valgrad = parallel::COSINEvalgrad;
	test[5].StartingGess = parallel::COSINEStartingGess;
	//7-CRAGGLVY
	//8-CURLY10
	//9-CURLY20
	//10-CURLY30
	//11-DIXMAANA
	//12-DIXMAANB
	//13-DIXMAANC
	//14-DIXMAAND
	//15-DIXMAANE
	//16-DIXMAANF
	//17-DIXMAANG
	//18-DIXMAANH
	//19-DIXMAANI
	//20-DIXMAANJ
	//21-DIXMAANK
	//22-DIXMAANL
	//23-DIXON3DQ
	test[22].size = 5000;
	strcpy(test[22].name, "DIXON3DQ");
	test[22].value = parallel::DIXON3DQvalue;
	test[22].grad = parallel::DIXON3DQgrad;
	test[22].valgrad = parallel::DIXON3DQvalgrad;
	test[22].StartingGess = parallel::DIXON3DQStartingGess;
	//24-DQDRTIC
	test[23].size = 10000;
	strcpy(test[23].name, "DQDRTIC");
	test[23].value = parallel::DQDRTICvalue;
	test[23].grad = parallel::DQDRTICgrad;
	test[23].valgrad = parallel::DQDRTICvalgrad;
	test[23].StartingGess = parallel::DQDRTICStartingGess;
	//25-DQRTIC
	test[24].size = 100000;
	strcpy(test[24].name, "DQRTIC");
	test[24].value = parallel::DQRTICvalue;
	test[24].grad = parallel::DQRTICgrad;
	test[24].valgrad = parallel::DQRTICvalgrad;
	test[24].StartingGess = parallel::DQRTICStartingGess;
	//26-EDENSCH
	test[25].size = 2000;
	strcpy(test[25].name, "EDENSCH");
	test[25].value = parallel::EDENSCHvalue;
	test[25].grad = parallel::EDENSCHgrad;
	test[25].valgrad = parallel::EDENSCHvalgrad;
	test[25].StartingGess = parallel::EDENSCHStartingGess;
	//27-EG2
	test[26].size = 5000;
	strcpy(test[26].name, "EG2");
	test[26].value = parallel::EG2value;
	test[26].grad = parallel::EG2grad;
	test[26].valgrad = parallel::EG2valgrad;
	test[26].StartingGess = parallel::EG2StartingGess;
	//28-ENGVAL1
	test[27].size = 5000;
	strcpy(test[27].name, "ENGVAL1");
	test[27].value = parallel::ENGVAL1value;
	test[27].grad = parallel::ENGVAL1grad;
	test[27].valgrad = parallel::ENGVAL1valgrad;
	test[27].StartingGess = parallel::ENGVAL1StartingGess;
	//EXTROSNB
	//FLETCHR
	test[29].size = 1000;
	strcpy(test[29].name, "FLETCHR");
	test[29].value = parallel::FLETCHRvalue;
	test[29].grad = parallel::FLETCHRgrad;
	test[29].valgrad = parallel::FLETCHRvalgrad;
	test[29].StartingGess = parallel::FLETCHRStartingGess;
	//31-FREUROTH
	test[30].size = 5000;
	strcpy(test[30].name, "FREUROTH");
	test[30].value = parallel::FREUROTHvalue;
	test[30].grad = parallel::FREUROTHgrad;
	test[30].valgrad = parallel::FREUROTHvalgrad;
	test[30].StartingGess = parallel::FREUROTHStartingGess;
	//32-GENHUMPS
	test[31].size = 5000;
	strcpy(test[31].name, "GENHUMPS");
	test[31].value = parallel::GENHUMPSvalue;
	test[31].grad = parallel::GENHUMPSgrad;
	test[31].valgrad = parallel::GENHUMPSvalgrad;
	test[31].StartingGess = parallel::GENHUMPSStartingGess;
	//33-GENROSE
	test[32].size = 1000;
	strcpy(test[32].name, "GENROSE");
	test[32].value = parallel::GENROSEvalue;
	test[32].grad = parallel::GENROSEgrad;
	test[32].valgrad = parallel::GENROSEvalgrad;
	test[32].StartingGess = parallel::GENROSEStartingGess;
	//34-LIARWDH
	test[33].size = 5000;
	strcpy(test[33].name, "LIARWDH");
	test[33].value = parallel::LIARWDHvalue;
	test[33].grad = parallel::LIARWDHgrad;
	test[33].valgrad = parallel::LIARWDHvalgrad;
	test[33].StartingGess = parallel::LIARWDHStartingGess;
	//35-MOREBV
	//36-NONCVXU2
	//37-NONDIA
	//38-NONDQUAR
	test[37].size = 10000;
	strcpy(test[37].name, "NONDQUAR");
	test[37].value = parallel::NONDQUARvalue;
	test[37].grad = parallel::NONDQUARgrad;
	test[37].valgrad = parallel::NONDQUARvalgrad;
	test[37].StartingGess = parallel::NONDQUARStartingGess;
	//39-PENALTY1
	//40-PENALTY2
	//41-POWER
	test[40].size = 1000;
	strcpy(test[40].name, "POWER");
	test[40].value = parallel::POWERvalue;
	test[40].grad = parallel::POWERgrad;
	test[40].valgrad = parallel::POWERvalgrad;
	test[40].StartingGess = parallel::POWERStartingGess;
	//42-QARTC
	//43-SROSENBR
	test[42].size = 10000;
	strcpy(test[42].name, "SROSENBR");
	test[42].value = parallel::SROSENBRvalue;
	test[42].grad = parallel::SROSENBRgrad;
	test[42].valgrad = parallel::SROSENBRvalgrad;
	test[42].StartingGess = parallel::SROSENBRStartingGess;
	//44-TRIDIA
	//45-WOODS
}
