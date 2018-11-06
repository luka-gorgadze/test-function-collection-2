#pragma once
INT n;			//The size of the problem

				/*===== 0 =================== ARWHEAD Function ================ 5000 ===================*/
double ARWHEADvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double group1 = 0.0;
	g[n - 1] = 0;
	for (i = 0; i < n - 1; i++)
	{
		group1 = pow(x[i], 2) + pow(x[n - 1], 2);
		fx += group1*group1 + 3.0 - 4.0*x[i];
		g[i] = -4.0 + 4.0*group1*x[i];
		g[n - 1] += 4.0*group1*x[n - 1];
	}
	return fx;
}
void ARWHEADgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double group1 = 0.0;
	g[n - 1] = 0;
	for (i = 0; i < n - 1; i++)
	{
		group1 = pow(x[i], 2) + pow(x[n - 1], 2);
		g[i] = -4.0 + 4.0*group1*x[i];
		g[n - 1] += 4.0*group1*x[n - 1];
	}
}
double ARWHEADvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double group1 = 0.0;
	for (i = 0; i < n - 1; i++)
	{
		group1 = pow(x[i], 2) + pow(x[n - 1], 2);
		fx += group1*group1 + 3.0 - 4.0*x[i];
	}
	return fx;
}
void ARWHEADStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 1.0;
}
/*==== 1 ================= BDQRTIC Function ============= 5000  ==================*/
double BDQRTICvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double group1 = 0.0, group2 = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 4; i++)
	{
		group1 = 3.0 - 4.0*x[i];
		group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
			+ 4 * pow(x[i + 3], 2) + 5 * pow(x[n - 1], 2);
		fx += group1*group1 + group2*group2;
		g[i] += -8.0*group1 + 4.0*group2*x[i];
		g[i + 1] += 8.0*group2*x[i + 1];
		g[i + 2] += 12.0*group2*x[i + 2];
		g[i + 3] += 16.0*group2*x[i + 3];
		g[n - 1] += 20.0*group2*x[n - 1];
	}
	return fx;
}
void BDQRTICgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double group1 = 0.0, group2 = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 4; i++)
	{
		group1 = 3.0 - 4.0*x[i];
		group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
			+ 4 * pow(x[i + 3], 2) + 5 * pow(x[n - 1], 2);
		g[i] += -8.0*group1 + 4.0*group2*x[i];
		g[i + 1] += 8.0*group2*x[i + 1];
		g[i + 2] += 12.0*group2*x[i + 2];
		g[i + 3] += 16.0*group2*x[i + 3];
		g[n - 1] += 20.0*group2*x[n - 1];
	}
}
double BDQRTICvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double group1 = 0.0, group2 = 0.0;

	for (i = 0; i < n - 4; i++)
	{
		group1 = 3.0 - 4.0*x[i];
		group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
			+ 4 * pow(x[i + 3], 2) + 5 * pow(x[n - 1], 2);
		fx += group1*group1 + group2*group2;
	}
	return fx;
}
void BDQRTICStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 1.0;
}
/*===== 2 ================= BROYDN7D Function ============= 5000 ==================*/
//precondition: n>=4
double BROYDN7Dvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double first = -2.0*x[1] + 1 + (3. - 2.0*x[0])*x[0];
	double last = -x[n - 2] + 1 + (3. - 2.0*x[n - 1])*x[n - 1];
	double fx = pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);
	for (i = 0; i < n; i++)
		g[i] = 0;

	g[0] = (7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1))*(3. - 4.*x[0]);
	g[1] = -(14.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
	g[n - 2] = -(7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
	g[n - 1] = (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1))*(3. - 4.*x[n - 1]);

	last = x[0] + x[n / 2];
	fx += pow(fabs(last), 7 / 3.0);
	g[0] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
	g[n / 2] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));

	for (i = 1; i < n / 2; i++)
	{
		first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];
		last = x[i] + x[i + n / 2];
		fx += pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);

		g[i - 1] += -(7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
		g[i] += (7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1))*(3 - 4 * x[i])
			+ (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
		g[i + 1] += -(14.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
		g[i + n / 2] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
	}
	for (i = n / 2; i < n - 1; i++)
	{
		first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];
		fx += pow(fabs(first), 7 / 3.0);

		g[i - 1] += -(7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
		g[i] += (7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1))*(3 - 4 * x[i]);
		g[i + 1] += -(14.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
	}
	return fx;
}
void BROYDN7Dgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double first = -2.0*x[1] + 1 + (3. - 2.0*x[0])*x[0];
	double last = -x[n - 2] + 1 + (3. - 2.0*x[n - 1])*x[n - 1];
	double fx = pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);
	for (i = 0; i < n; i++)
		g[i] = 0;

	g[0] = (7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1))*(3. - 4.*x[0]);
	g[1] = -(14.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
	g[n - 2] = -(7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
	g[n - 1] = (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1))*(3. - 4.*x[n - 1]);

	last = x[0] + x[n / 2];
	g[0] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
	g[n / 2] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));

	for (i = 1; i < n / 2; i++)
	{
		first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];
		last = x[i] + x[i + n / 2];

		g[i - 1] += -(7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
		g[i] += (7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1))*(3 - 4 * x[i])
			+ (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
		g[i + 1] += -(14.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
		g[i + n / 2] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last>0) ? (1) : (-1));
	}
	for (i = n / 2; i < n - 1; i++)
	{
		first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];

		g[i - 1] += -(7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
		g[i] += (7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1))*(3 - 4 * x[i]);
		g[i + 1] += -(14.0 / 3)* pow(fabs(first), 4 / 3.0)*((first>0) ? (1) : (-1));
	}
}
double BROYDN7Dvalue
(
	double *x,
	INT n
)
{
	INT i;
	double first = -2.0*x[1] + 1 + (3. - 2.0*x[0])*x[0];
	double last = -x[n - 2] + 1 + (3. - 2.0*x[n - 1])*x[n - 1];
	double fx = pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);
	last = x[0] + x[n / 2];
	fx += pow(fabs(last), 7 / 3.0);

	for (i = 1; i < n / 2; i++)
	{
		first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];
		last = x[i] + x[i + n / 2];
		fx += pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);
	}
	for (i = n / 2; i < n - 1; i++)
	{
		first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];
		fx += pow(fabs(first), 7 / 3.0);
	}
	return fx;
}
void BROYDN7DStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 1.0;
}
/*===== 3 =================== BRYBND Function ================ 5000 ===================*/
double BRYBNDvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	double group = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	int mL = 5, mU = 1;
	for (i = 0; i < n - 1; i++)
	{
		double bnd = 0.0;
		int  lo = (0 > i - mL) ? (0) : (i - mL), u = (n - 1 < i + mU) ? (n - 1) : (i + mU);
		for (j = lo; j < i; j++)
			bnd += x[j] * (1 + x[j]);
		for (j = i + 1; j <= u; j++)
			bnd += x[j] * (1 + x[j]);
		group = x[i] * (2.0 + 5.0*pow(x[i], 2.0)) + 1 - bnd;
		fx += pow(group, 2.0);
		g[i] += 2.0*group*(2.0 + 15.0*pow(x[i], 2.0));
		for (j = lo; j < i; j++)
			g[j] -= 2.0*group*(1 + 2.0*x[j]);
		for (j = i + 1; j <= u; j++)
			g[j] -= 2.0*group*(1 + 2.0*x[j]);
	}
	return fx;
}
void BRYBNDgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	double group = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	int mL = 5, mU = 1;
	for (i = 0; i < n - 1; i++)
	{
		double bnd = 0.0;
		int  lo = (0 > i - mL) ? (0) : (i - mL), u = (n - 1 < i + mU) ? (n - 1) : (i + mU);
		for (j = lo; j < i; j++)
			bnd += x[j] * (1 + x[j]);
		for (j = i + 1; j <= u; j++)
			bnd += x[j] * (1 + x[j]);
		group = x[i] * (2.0 + 5.0*pow(x[i], 2.0)) + 1 - bnd;
		g[i] += 2.0*group*(2.0 + 15.0*pow(x[i], 2.0));
		for (j = lo; j < i; j++)
			g[j] -= 2.0*group*(1 + 2.0*x[j]);
		for (j = i + 1; j <= u; j++)
			g[j] -= 2.0*group*(1 + 2.0*x[j]);
	}
}
double BRYBNDvalue
(
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	double group = 0.0;
	int mL = 5, mU = 1;
	for (i = 0; i < n - 1; i++)
	{
		double bnd = 0.0;
		int  lo = (0 > i - mL) ? (0) : (i - mL), u = (n - 1 < i + mU) ? (n - 1) : (i + mU);
		for (j = lo; j < i; j++)
			bnd += x[j] * (1 + x[j]);
		for (j = i + 1; j <= u; j++)
			bnd += x[j] * (1 + x[j]);
		group = x[i] * (2.0 + 5.0*pow(x[i], 2.0)) + 1 - bnd;
		fx += pow(group, 2.0);
	}
	return fx;
}
void BRYBNDStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = -1.0;
}
/*========= 4 ========= Chained Wood - CHAINWOO ======== 4000 ==========*/
double CHAINWOOvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1, item2, item3, item4, item5, item6;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 3; i += 2)
	{
		item1 = (x[i + 1] - pow(x[i], 2.0));
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - pow(x[i + 2], 2.0));
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		fx += 100 * pow(item1, 2.0) + pow(item2, 2.0) + 90 * pow(item3, 2.0)
			+ pow(item4, 2.0) + 10.0*pow(item5, 2.0) + 0.1*pow(item6, 2.0);
		g[i] -= (400 * item1*x[i] + 2 * item2);
		g[i + 1] += 200 * item1 + 20.0*item5 + 0.2*item6;
		g[i + 2] -= (360 * item3*x[i + 2] + 2.0*item4);
		g[i + 3] += 180 * item3 + 20.0*item5 - 0.2*item6;
	}
	return fx;
}
void CHAINWOOgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1, item2, item3, item4, item5, item6;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 3; i += 2)
	{
		item1 = (x[i + 1] - pow(x[i], 2.0));
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - pow(x[i + 2], 2.0));
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		g[i] -= (400 * item1*x[i] + 2 * item2);
		g[i + 1] += 200 * item1 + 20.0*item5 + 0.2*item6;
		g[i + 2] -= (360 * item3*x[i + 2] + 2.0*item4);
		g[i + 3] += 180 * item3 + 20.0*item5 - 0.2*item6;
	}
}
double CHAINWOOvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1, item2, item3, item4, item5, item6;

	for (i = 0; i < n - 3; i += 2)
	{
		item1 = (x[i + 1] - pow(x[i], 2.0));
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - pow(x[i + 2], 2.0));
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		fx += 100 * pow(item1, 2.0) + pow(item2, 2.0) + 90 * pow(item3, 2.0)
			+ pow(item4, 2.0) + 10.0*pow(item5, 2.0) + 0.1*pow(item6, 2.0);
	}
	return fx;
}
void CHAINWOOStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	x[0] = -3.0;
	x[1] = -1.0;
	x[2] = -3.0;
	x[3] = -1.0;
	for (i = 4; i<n; i++)
		x[i] = -2.0;
}
/*========= 5 ========== COSINE Function =========== 10000 ===========*/
double COSINEvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = 0.0;

	for (i = 0; i < n; i++)
		g[i] = 0.0;

	for (i = 0; i < n - 1; i++)
	{
		item = -0.5*x[i + 1] + x[i] * x[i];
		fx += cos(item);
		g[i] -= 2.0*sin(item)*x[i];
		g[i + 1] += 0.5*sin(item);
	}
	return fx;
}
void COSINEgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;

	for (i = 0; i < n; i++)
		g[i] = 0.0;

	for (i = 0; i < n - 1; i++)
	{
		item = -0.5*x[i + 1] + x[i] * x[i];
		g[i] -= 2.0*sin(item)*x[i];
		g[i + 1] += 0.5*sin(item);
	}
}
double COSINEvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = 0.0;

	for (i = 0; i < n - 1; i++)
	{
		item = -0.5*x[i + 1] + x[i] * x[i];
		fx += cos(item);
	}
	return fx;
}
void COSINEStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i++)
		x[i] = 1.0;
}
/*============ 6 ========= CRAGGLVY  =========== 5000 =============*/
double CRAGGLVYvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1, item2, element, item3, item4;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 3; i += 2)
	{
		item1 = exp(x[i]) - x[i + 1];
		item2 = x[i + 1] - x[i + 2];
		element = x[i + 2] - x[i + 3];
		item3 = tan(element) + element;
		item4 = (x[i + 3] - 1);

		fx += pow(item1, 4.0) + 100 * pow(item2, 6.0)
			+ pow(item3, 4.0) + pow(x[i], 8.0) + pow(item4, 2.0);
		g[i] += 4 * pow(item1, 3.0)*exp(x[i]) + 8 * pow(x[i], 7.0);
		g[i + 1] += -4 * pow(item1, 3.0) + 600 * pow(item2, 5.0);
		g[i + 2] += -600 * pow(item2, 5.0)
			+ 4 * pow(item3, 3.0)*(1 / pow(cos(element), 2.0) + 1);

		g[i + 3] += -4 * pow(item3, 3.0)*(1 / pow(cos(element), 2.0) + 1)
			+ 2 * item4;
	}
	return fx;
}
void CRAGGLVYgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2, element, item3, item4;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 3; i += 2)
	{
		item1 = exp(x[i]) - x[i + 1];
		item2 = x[i + 1] - x[i + 2];
		element = x[i + 2] - x[i + 3];
		item3 = tan(element) + element;
		item4 = (x[i + 3] - 1);


		g[i] += 4 * pow(item1, 3.0)*exp(x[i]) + 8 * pow(x[i], 7.0);
		g[i + 1] += -4 * pow(item1, 3.0) + 600 * pow(item2, 5.0);
		g[i + 2] += -600 * pow(item2, 5.0)
			+ 4 * pow(item3, 3.0)*(1 / pow(cos(element), 2.0) + 1);

		g[i + 3] += -4 * pow(item3, 3.0)*(1 / pow(cos(element), 2.0) + 1)
			+ 2 * item4;
	}
}
double CRAGGLVYvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1, item2, element, item3, item4;

	for (i = 0; i < n - 3; i += 2)
	{
		item1 = exp(x[i]) - x[i + 1];
		item2 = x[i + 1] - x[i + 2];
		element = x[i + 2] - x[i + 3];
		item3 = tan(element) + element;
		item4 = (x[i + 3] - 1);

		fx += pow(item1, 4.0) + 100 * pow(item2, 6.0)
			+ pow(item3, 4.0) + pow(x[i], 8.0) + pow(item4, 2.0);
	}
	return fx;
}
void CRAGGLVYStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	x[0] = 1.0;

	for (i = 1; i<n; i++)
		x[i] = 2.0;
}
/*===== 7 =============== CURLY10 Function ============= 500 ===================*/
double CURLY10valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	int k = 10;
	int ipk;		// min{i+k,n-1}
	double q;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
		for (j = i; j <= ipk; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;
	}
	return fx;
}
void CURLY10grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	int k = 10;
	int ipk;		// min{i+k,n-1}
	double q;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		for (j = i; j <= ipk; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;
	}
}
double CURLY10value
(
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	int k = 10;
	int ipk;		// min{i+k,n-1}
	double q;

	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
	}
	return fx;
}
void CURLY10StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 0.0001 / (n + 1);
}
/*===== 8 =============== CURLY20 Function ============= 500 ===================*/
double CURLY20valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	int k = 20;
	int ipk;		// min{i+k,n-1}
	double q;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
		for (j = i; j <= ipk; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;
	}
	return fx;
}
void CURLY20grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	int k = 20;
	int ipk;		// min{i+k,n-1}
	double q;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		for (j = i; j <= ipk; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;
	}
}
double CURLY20value
(
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	int k = 20;
	int ipk;		// min{i+k,n-1}
	double q;

	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
	}
	return fx;
}
void CURLY20StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 0.0001 / (n + 1);
}
/*===== 9 =============== CURLY10 Function ============= 500 ===================*/
double CURLY30valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	int k = 30;
	int ipk;		// min{i+k,n-1}
	double q;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
		for (j = i; j <= ipk; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;
	}
	return fx;
}
void CURLY30grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	int k = 30;
	int ipk;		// min{i+k,n-1}
	double q;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		for (j = i; j <= ipk; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;
	}
}
double CURLY30value
(
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	int k = 30;
	int ipk;		// min{i+k,n-1}
	double q;

	for (i = 0; i < n; i++)
	{
		q = 0.0;
		ipk = (((i + k) < (n - 1)) ? (i + k) : (n - 1));
		for (j = i; j <= ipk; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
	}
	return fx;
}
void CURLY30StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 0.0001 / (n + 1);
}
/*===== 10 ======================= DIXMAANA Function ================== 5000 ================*/
double DIXMAANAvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*pow(item, 2) + 0.125*x[i] * x[i + 2 * m];
		g[i] += 2 * x[i] + 0.25*item *pow(x[i + m], 2) + 0.125*x[i + 2 * m];
		g[i + m] += 0.5*item *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*pow(item, 2);
		g[i] += 2 * x[i] + 0.25*item *pow(x[i + m], 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n; i++)
	{
		fx += pow(x[i], 2);
		g[i] += 2 * x[i];
	}
	return fx;
}
void DIXMAANAgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.25*item *pow(x[i + m], 2) + 0.125*x[i + 2 * m];
		g[i + m] += 0.5*item *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.25*item *pow(x[i + m], 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n; i++)
	{
		g[i] += 2 * x[i];
	}
}
double DIXMAANAvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item = 0.0;
	int m = n / 3;

	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*pow(item, 2) + 0.125*x[i] * x[i + 2 * m];
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*pow(item, 2);
	}
	for (i = 2 * m; i < n; i++)
	{
		fx += pow(x[i], 2);
	}
	return fx;
}
void DIXMAANStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 2.0;
}
/*=====11  ======================= DIXMAANB Function ================== 5000 ================*/
double DIXMAANBvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + item1*item1 + item2*item2 + (0.0625*x[i] * x[i + 2 * m]);
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2) +
			(0.0625*x[i + 2 * m]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.0625*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + item1*item1 + item2*item2;
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + item1*item1;
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANBgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2) +
			(0.0625*x[i + 2 * m]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.0625*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANBvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;

	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + item1*item1 + item2*item2 + (0.0625*x[i] * x[i + 2 * m]);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + item1*item1 + item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*=====12  ======================= DIXMAANC Function ================== 5000 ================*/
double DIXMAANCvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*item1*item1 + 0.125*item2*item2 + 0.125*x[i] * x[i + 2 * m];
		g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2) + 0.125*x[i + 2 * m];
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*item1*item1 + 0.125*item2*item2;
		g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + 0.125*item1*item1;
		g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANCgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2) + 0.125*x[i + 2 * m];
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANCvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;

	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*item1*item1 + 0.125*item2*item2 + 0.125*x[i] * x[i + 2 * m];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.125*item1*item1 + 0.125*item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + 0.125*item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*===== 13 ======================= DIXMAAND Function ================== 5000 ================*/
double DIXMAANDvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.26*item1*item1 + 0.26*item2*item2 + 0.26*x[i] * x[i + 2 * m];
		g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2) + 0.26*x[i + 2 * m];
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.26*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.26*item1*item1 + 0.26*item2*item2;
		g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + 0.26*item1*item1;
		g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void  DIXMAANDgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2) + 0.26*x[i + 2 * m];
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.26*x[i];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANDvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;

	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.26*item1*item1 + 0.26*item2*item2 + 0.26*x[i] * x[i + 2 * m];
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2) + 0.26*item1*item1 + 0.26*item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2) + 0.26*item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*===== 14 ======================= DIXMAANE Function ================== 3000 ================*/
double DIXMAANEvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*pow(item, 2) + (0.125*x[i] * x[i + 2 * m] * (i + 1)) / n;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item *pow(x[i + m], 2) + (0.125*x[i + 2 * m] * (i + 1)) / n;
		g[i + m] += 0.5*item *x[i] * x[i + m];
		g[i + 2 * m] += (0.125*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*pow(item, 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item *pow(x[i + m], 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n; i++)
	{
		fx += (pow(x[i], 2)*i) / n;
		g[i] += (2 * x[i] * (i + 1)) / n;
	}
	return fx;
}
void DIXMAANEgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item *pow(x[i + m], 2) + (0.125*x[i + 2 * m] * (i + 1)) / n;
		g[i + m] += 0.5*item *x[i] * x[i + m];
		g[i + 2 * m] += (0.125*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item *pow(x[i + m], 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n; i++)
	{
		g[i] += (2 * x[i] * (i + 1)) / n;
	}
}
double DIXMAANEvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item = 0.0;
	int m = n / 3;
	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*pow(item, 2) + (0.125*x[i] * x[i + 2 * m] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*pow(item, 2);
	}
	for (i = 2 * m; i < n; i++)
		fx += (pow(x[i], 2)*i) / n;
	return fx;
}
/*===== 15 ======================= DIXMAANF Function ================== 3000 ================*/
double DIXMAANFvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + item1*item1 + item2*item2 + (0.0625*x[i] * x[i + 2 * m] * (i + 1)) / n;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2) +
			(0.0625*x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
		g[i + 2 * m] += (0.0625*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + item1*item1 + item2*item2;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2)*(i + 1)) / n + item1*item1;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANFgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2) +
			(0.0625*x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
		g[i + 2 * m] += (0.0625*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANFvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = n / 3;
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + item1*item1 + item2*item2 + (0.0625*x[i] * x[i + 2 * m] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + item1*item1 + item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2)*(i + 1)) / n + item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*===== 16 ======================= DIXMAANG Function ================== 3000 ================*/
double DIXMAANGvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1 + 0.125*item2*item2
			+ (0.125*x[i] * x[i + 2 * m] * (i + 1)) / n;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2)
			+ (0.125*x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
		g[i + 2 * m] += (0.125*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1 + 0.125*item2*item2;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANGgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2)
			+ (0.125*x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
		g[i + 2 * m] += (0.125*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANGvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = (n / 3);

	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1 + 0.125*item2*item2
			+ (0.125*x[i] * x[i + 2 * m] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1 + 0.125*item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*===== 17 ======================= DIXMAANH Function ================== 3000 ================*/
double DIXMAANHvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1 + 0.26*item2*item2
			+ (0.26*x[i] * x[i + 2 * m] * (i + 1)) / n;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2)
			+ (0.26*x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
		g[i + 2 * m] += (0.26*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1 + 0.26*item2*item2;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1;
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANHgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2)
			+ (0.26*x[i + 2 * m] * (i + 1)) / n;
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
		g[i + 2 * m] += (0.26*x[i] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANHvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = 0.0;
	double item2 = 0.0;
	int m = (n / 3);

	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1 + 0.26*item2*item2
			+ (0.26*x[i] * x[i + 2 * m] * (i + 1)) / n;
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1 + 0.26*item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*===== 18 ======================= DIXMAANI Function ================== 3000 ================*/
double DIXMAANIvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item = 0.0;
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*pow(item, 2)
			+ 0.125*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item *pow(x[i + m], 2)
			+ 0.125*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*pow(item, 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item *pow(x[i + m], 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n; i++)
	{
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2);
	}
	return fx;
}
void DIXMAANIgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = 0.0;
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item *pow(x[i + m], 2)
			+ 0.125*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item *pow(x[i + m], 2);
		g[i + m] += 0.5*item *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n; i++)
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2);
}
double DIXMAANIvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item = 0.0;
	int m = (n / 3);

	for (i = 0; i < m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*pow(item, 2)
			+ 0.125*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*pow(item, 2);
	}
	for (i = 2 * m; i < n; i++)
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2);

	return fx;
}
/*===== 19 ======================= DIXMAANJ Function ================== 3000 ================*/
double DIXMAANJvalgrad
(
	double *g,
	double *x,
	INT nn
)
{
	INT i;
	double fx = 1.0;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + item1*item1 + item2*item2
			+ 0.0625*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5*item2*pow(x[i + m], 2) + 0.0625*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.0625*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + item1*item1 + item2*item2;
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + item1*item1;
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANJgrad
(
	double *g,
	double *x,
	INT nn
)
{
	INT i;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5*item2*pow(x[i + m], 2) + 0.0625*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.0625*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANJvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + item1*item1 + item2*item2
			+ 0.0625*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = 0.25* x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + item1*item1 + item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*===== 20 ======================= DIXMAANK Function ================== 5000 ================*/
double DIXMAANKvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*item1*item1 + 0.125*item2*item2
			+ 0.125*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.25*item2*pow(x[i + m], 2) + 0.125*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*item1*item1 + 0.125*item2*item2;
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.25*item2*pow(x[i + m], 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*item1*item1;
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANKgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.25*item2*pow(x[i + m], 2) + 0.125*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.125*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.25*item2*pow(x[i + m], 2);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 0.5* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANKvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);

	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*item1*item1 + 0.125*item2*item2
			+ 0.125*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*item1*item1 + 0.125*item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.125*item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*===== 21 ======================= DIXMAANL Function ================== 5000 ================*/
double DIXMAANLvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.26*item1*item1 + 0.26*item2*item2
			+ 0.26*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.52*item2*pow(x[i + m], 2) + 0.26*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.26*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.26*item1*item1 + 0.26*item2*item2;
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.52*item2*pow(x[i + m], 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.26*item1*item1;
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	fx += pow(x[n - 1], 2);
	g[n - 1] += 2 * x[n - 1];
	return fx;
}
void DIXMAANLgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.52*item2*pow(x[i + m], 2) + 0.26*x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
		g[i + 2 * m] += 0.26*x[i] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.52*item2*pow(x[i + m], 2);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += 1.04* item2 *x[i] * x[i + m];
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] * pow((i + 1) / (double)(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
	}
	g[n - 1] += 2 * x[n - 1];
}
double DIXMAANLvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1 = (0.0);
	double item2 = (0.0);
	int m = (n / 3);
	for (i = 0; i < m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.26*item1*item1 + 0.26*item2*item2
			+ 0.26*x[i] * x[i + 2 * m] * pow((i + 1) / (double)(n), 2);
	}
	for (i = m; i < 2 * m; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		item2 = x[i] * pow(x[i + m], 2);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.26*item1*item1 + 0.26*item2*item2;
	}
	for (i = 2 * m; i < n - 1; i++)
	{
		item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		fx += pow(x[i], 2)*pow((i + 1) / (double)(n), 2) + 0.26*item1*item1;
	}
	fx += pow(x[n - 1], 2);
	return fx;
}
/*=============== 22 ============= DIXON3DQ Function ============= 5000 =================*/
double DIXON3DQvalgrad
(
	double *g,
	double *x,
	INT n
)
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
void DIXON3DQgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	item = x[0] - 1;
	g[0] = 2.0*item;

	for (i = 1; i < n - 1; i++)
		g[i] = 0.0;

	item = x[n - 1] - 1;
	g[n - 1] = 2.0*item;

	for (i = 1; i < n - 1; i++)
	{
		item = (x[i] - x[i + 1]);
		g[i] += 2.0*item;
		g[i + 1] -= 2.0*item;
	}
}
double DIXON3DQvalue
(
	double *x,
	INT n
)
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
void DIXON3DQStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i++)
		x[i] = -1.0;
}
/*===== 23 ===================== DQDRTIC Function ============= 10000 =================*/
double DQDRTICvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 2; i++)
	{
		fx += pow(x[i], 2) + 100 * pow(x[i + 1], 2) + 100 * pow(x[i + 2], 2);
		g[i] += 2 * x[i];
		g[i + 1] += 200 * x[i + 1];
		g[i + 2] += 200 * x[i + 2];
	}
	return fx;
}
void DQDRTICgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 2; i++)
	{
		g[i] += 2 * x[i];
		g[i + 1] += 200 * x[i + 1];
		g[i + 2] += 200 * x[i + 2];
	}
}
double DQDRTICvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	for (i = 0; i < n - 2; i++)
	{
		fx += pow(x[i], 2) + 100 * pow(x[i + 1], 2) + 100 * pow(x[i + 2], 2);
	}
	return fx;
}
void DQDRTICStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i <n; i++)
		x[i] = 3.0;
}
/*===== 24 ===================== DQRTIC Function ===== 5000 ======================*/
/*
//valgrad:
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
double DQRTICvalgradRangedParallel
(
	double *g,
	double *x,
	INT begin,
	INT end
)
{
	if ((end - begin)  < min_per_thread)
		return DQRTICvalgradRanged(g, x, begin, end);
	std::future<double> first_half_result =
		std::async(DQRTICvalgradRangedParallel, g, x, begin,(begin+end)/ 2);
	double second_half_result = DQRTICvalgradRangedParallel(g, x, (begin + end) / 2,end);
	return first_half_result.get() + second_half_result;
}
double DQRTICvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	return DQRTICvalgradRangedParallel(g, x, 0, n);
}
//grad:
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
void DQRTICgradRangedParallel
(
	double *g,
	double *x,
	INT begin,
	INT end
)
{
	if ((end - begin)  < min_per_thread)
		return DQRTICgradRanged(g, x, begin, end);
	auto first_half_result = std::async(DQRTICgradRangedParallel, g, x, begin, (begin + end) / 2);
	DQRTICgradRangedParallel(g, x, (begin + end) / 2, end);
	first_half_result.get();
}


void DQRTICgrad
(
	double *g,
	double *x,
	INT n
)
{
	DQRTICgradRangedParallel(g, x, 0, n);
}
//value:
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
double DQRTICvalueRangedParallel
(
	double *x,
	INT begin,
	INT end
)
{
	if ((end - begin)  < min_per_thread)
		return DQRTICvalueRanged(x, begin, end);
	std::future<double> first_half_result =
		std::async(DQRTICvalueRangedParallel, x, begin, (begin + end) / 2);
	double second_half_result = DQRTICvalueRangedParallel(x, (begin + end) / 2, end);
	return first_half_result.get() + second_half_result;
}
double DQRTICvalue
(
	double *x,
	INT n
)
{
	return DQRTICvalueRangedParallel(x, 0, n);
}
void DQRTICStartingGessRanged
(
	double *x,
	INT begin,
	INT end
)
{
	INT i;
	for (i = begin; i < end; i++)
	{
		x[i] = 2.0;
	}
}
void DQRTICStartingGessRangedParallel
(
	double *x,
	INT begin,
	INT end
)
{
	if ((end- begin) < min_per_thread)
		return DQRTICStartingGessRanged(x, begin, end);
	auto first_half_result = std::async(DQRTICStartingGessRangedParallel, x, begin, (begin + end) / 2);
	DQRTICStartingGessRangedParallel(x, (begin + end) / 2, end);
	first_half_result.get();
}
void DQRTICStartingGess
(
	double *x,
	INT n
)
{
	 DQRTICStartingGessRangedParallel(x, 0, n);
}
*/
/*

void DQRTICvalgradBlock
(
	double *g,
	double *x,
	INT begin,
	INT end,
	double *res
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
	*res = f;
}
double DQRTICvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT block_start = 0,i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		threads[i] = std::thread(DQRTICvalgradBlock, g, x, block_start, block_end, &results[i]);
		block_start = block_end;
	}
	DQRTICvalgradBlock(g, x, block_start, n, &results[num_threads-1]);
	for (auto i = threads.begin(); i != threads.end(); ++i) i->join();
	return std::accumulate(results.begin(), results.end(), 0);
}

void DQRTICgradBlock
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
	INT block_start = 0,i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		threads[i] = std::thread(DQRTICgradBlock, g, x, block_start, block_end);
		block_start = block_end;
	}
	DQRTICgradBlock(g, x, block_start, n);
	for (auto i = threads.begin(); i != threads.end(); ++i) i->join();
}

void DQRTICvalueBlock
(
	double *x,
	INT begin,
	INT end,
	double *res
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
	*res = f;
}
double DQRTICvalue
(
	double *x,
	INT n
)
{
	INT block_start = 0,i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		threads[i] = std::thread(DQRTICvalueBlock, x, block_start, block_end, &results[i]);
		block_start = block_end;
	}
	DQRTICvalueBlock(x, block_start, n, &results[num_threads - 1]);
	for (auto i = threads.begin(); i != threads.end(); ++i) i->join();
	return std::accumulate(results.begin(), results.end(), 0);
}
void DQRTICStartingGessBlock
(
	double *x,
	INT begin,
	INT end
)
{
	INT i;
	for (i = begin; i < end; i++)
	{
		x[i] = 2.0;
	}
}
void DQRTICStartingGess
(
	double *x,
	INT n
)
{
	INT block_start = 0, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		threads[i] = std::thread(DQRTICStartingGessBlock, x, block_start, block_end);
		block_start = block_end;
	}
	DQRTICStartingGessBlock(x, block_start, n);
	for (auto i = threads.begin(); i != threads.end(); ++i) i->join();
}
*/
//-----------------------------------------------------------------------
double DQRTICvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double item;
	for (i = 0; i < n; i++)
	{
		item = x[i] - i - 1;
		fx += pow(item, 4);
		g[i] = 4 * pow(item, 3);
	}
	return fx;
}
void DQRTICgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	for (i = 0; i < n; i++)
	{
		item = x[i] - i - 1;
		g[i] = 4 * pow(item, 3);
	}
}
double DQRTICvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double item;
	for (i = 0; i < n; i++)
	{
		item = x[i] - i - 1;
		fx += pow(item, 4);
	}
	return fx;
}
void DQRTICStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i <n; i++)
		x[i] = 2.0;
}


/*===== 25 ===================== EDENSCH Function ===== 2000 ======================*/
double EDENSCHvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2, item3;
	for (i = 0; i < n - 1; i++)
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
void EDENSCHgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2, item3;
	for (i = 0; i < n - 1; i++)
	{
		item1 = x[i] - 2;
		item2 = x[i] * x[i + 1] - 2 * x[i + 1];
		item3 = x[i + 1] + 1;
		g[i] += 4 * pow(item1, 3) + 2 * item2*x[i + 1];
		g[i + 1] += 2 * item2*(x[i] - 2.0) + 2 * item3;
	}
}
double EDENSCHvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double item1, item2, item3;
	for (i = 0; i < n - 1; i++)
	{
		item1 = x[i] - 2;
		item2 = x[i] * x[i + 1] - 2 * x[i + 1];
		item3 = x[i + 1] + 1;
		fx += 16 + pow(item1, 4) + pow(item2, 2) + pow(item3, 2);
	}
	return fx;
}
void EDENSCHStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i <n; i++)
		x[i] = 0.0;
}
/*===== 26 ===================== EG2 Function =====5000 ======================*/
double EG2valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		g[i] = 0;
	double fx = 0.5*sin(pow(x[n - 1], 2));
	g[n - 1] = cos(pow(x[n - 1], 2))*x[n - 1];
	double item;
	for (i = 0; i < n - 1; i++)
	{
		item = x[0] + pow(x[i], 2) - 1;;
		fx += sin(item);
		g[0] += cos(item);
		g[i] += 2 * cos(item)*x[i];
	}
	return fx;
}
void EG2grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		g[i] = 0;
	g[n - 1] = cos(pow(x[n - 1], 2))*x[n - 1];
	double item;
	for (i = 0; i < n - 1; i++)
	{
		item = x[0] + pow(x[i], 2) - 1;;
		g[0] += cos(item);
		g[i] += 2 * cos(item)*x[i];
	}
}
double EG2value
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.5*sin(pow(x[n - 1], 2));
	double item;
	for (i = 0; i < n - 1; i++)
	{
		item = x[0] + pow(x[i], 2) - 1;;
		fx += sin(item);
	}
	return fx;
}
void EG2StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i <n; i++)
		x[i] = 0.;
}
/*========= 27 ========== ENGVAL1 Function ============= 5000 ==============*/
double ENGVAL1valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = (0.0);

	for (i = 0; i < n; i++)
		g[i] = 0.0;

	for (i = 0; i < n - 1; i++)
	{
		item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
		fx += item*item + (3 - 4.0*x[i]);
		g[i] += 4.0*item*x[i] - 4.0;
		g[i + 1] += 4.0*item*x[i + 1];
	}
	return fx;
}
void ENGVAL1grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	for (i = 0; i < n; i++)
		g[i] = 0.0;

	for (i = 0; i < n - 1; i++)
	{
		item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
		g[i] += 4.0*item*x[i] - 4.0;
		g[i + 1] += 4.0*item*x[i + 1];
	}
}
double ENGVAL1value
(
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = (0.0);
	for (i = 0; i < n - 1; i++)
	{
		item = pow(x[i], 2.0) + pow(x[i + 1], 2.0);
		fx += item*item + (3 - 4.0*x[i]);
	}
	return fx;
}
void ENGVAL1StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i <n; i++)
		x[i] = 2.0;

}
/*===== 28 ======================= EXTROSNB ================== 1000 ================*/
double EXTROSNBvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1);
	double fx = (pow(item, 2));
	for (i = 1; i < n; i++)
		g[i] = 0.0;

	g[0] = 2 * item;
	for (i = 1; i < n; i++) {
		item = 10 * (x[i] - pow(x[i - 1], 2));
		fx += item * item;
		g[i - 1] += -40.0 * item *x[i - 1];
		g[i] += 20.0 * item;
	}
	return fx;
}
void EXTROSNBgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1);
	for (i = 1; i < n; i++)
		g[i] = 0.0;

	g[0] = 2 * item;
	for (i = 1; i < n; i++) {
		item = 10 * (x[i] - pow(x[i - 1], 2));
		g[i - 1] += -40.0 * item *x[i - 1];
		g[i] += 20.0 * item;
	}
}
double EXTROSNBvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1);
	double fx = (pow(item, 2));
	for (i = 1; i < n; i++) {
		item = 10 * (x[i] - pow(x[i - 1], 2));
		fx += item * item;
	}
	return fx;
}
void EXTROSNBStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i <n; i++)
		x[i] = -1.0;
}
/*========= 29 ========== FLETCHR Function ============= 1000 ==============*/
double FLETCHRvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;

	for (i = 0; i < n - 1; i++)
	{
		item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
		fx += item*item;
		g[i] += 20.0*item*(-2.0*x[i] - 1.0);
		g[i + 1] += 20.0*item;
	}
	return 100.*fx;
}
void FLETCHRgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;

	for (i = 0; i < n; i++)
		g[i] = 0;

	for (i = 0; i < n - 1; i++)
	{
		item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
		g[i] += 20.0*item*(-2.0*x[i] - 1.0);
		g[i + 1] += 20.0*item;
	}
}
double FLETCHRvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = 0.0;
	for (i = 0; i < n - 1; i++)
	{
		item = x[i + 1] - x[i] + 1 - pow(x[i], 2);
		fx += item*item;
	}
	return 100.*fx;
}
void FLETCHRStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i <n; i++)
		x[i] = 0.0;

}
/*========= 30 ==============  Freudenstein and Roth: FREUROTH ====== 5000 =============*/
double FREUROTHvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double item1, item2;;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 1; i++)
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
void FREUROTHgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 0; i < n - 1; i++)
	{
		item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
		item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
		g[i] += 2.0*item1 + 2.0*item2;
		g[i + 1] += 2.0*item1*(10 * x[i + 1] - 3.0*x[i + 1] * x[i + 1] - 2.0) +
			2.0*item2*(2 * x[i + 1] + 3.0*x[i + 1] * x[i + 1] - 14.0);
	}
}
double FREUROTHvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double item1, item2;;
	for (i = 0; i < n - 1; i++)
	{
		item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
		item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
		fx += item1*item1 + item2*item2;
	}
	return fx;
}
void FREUROTHStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	x[0] = 0.5;
	x[1] = -2.0;
	for (i = 2; i < n; i++)
		x[i] = 0.0;
}
/*========= 31 ===================  GENHUMPS ====== 5000 ======================*/
double GENHUMPSvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2;;
	for (i = 0; i < n - 1; i++)
	{

		item1 = sin(2.0*x[i]);
		item2 = sin(2.0*x[i + 1]);
		fx += item1*item1 *item2*item2 + 0.05*(pow(x[i], 2) + pow(x[i + 1], 2));
		g[i] += 4.0*item1*cos(2.0*x[i]) *item2*item2 + 0.1*x[i];
		g[i + 1] += 4.0*item2*cos(2.0*x[i + 1])* item1*item1 + 0.1*x[i + 1];
	}
	return fx;
}
void GENHUMPSgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2;;
	for (i = 0; i < n - 1; i++)
	{

		item1 = sin(2.0*x[i]);
		item2 = sin(2.0*x[i + 1]);
		g[i] += 4.0*item1*cos(2.0*x[i]) *item2*item2 + 0.1*x[i];
		g[i + 1] += 4.0*item2*cos(2.0*x[i + 1])* item1*item1 + 0.1*x[i + 1];
	}
}
double GENHUMPSvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	double item1, item2;;
	for (i = 0; i < n - 1; i++)
	{

		item1 = sin(2.0*x[i]);
		item2 = sin(2.0*x[i + 1]);
		fx += item1*item1 *item2*item2 + 0.05*(pow(x[i], 2) + pow(x[i + 1], 2));
	}
	return fx;
}
void GENHUMPSStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	x[0] = -506.0;
	for (i = 1; i < n; i++)
		x[i] = 506.2;
}
/*===== 32 ======================= GENROSE ====================== 1000 ===================*/
double GENROSEvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1, item2;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 1; i < n; i++)
	{
		item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
		item2 = x[i] - 1.0;
		fx += item1*item1 + item2*item2;
		g[i - 1] -= 40.0*item1 * x[i - 1];
		g[i] += 20 * item1 + 2 * item2;
	}
	return fx;
}
void GENROSEgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;
	for (i = 0; i < n; i++)
		g[i] = 0;
	for (i = 1; i < n; i++)
	{
		item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
		item2 = x[i] - 1.0;
		g[i - 1] -= 40.0*item1 * x[i - 1];
		g[i] += 20 * item1 + 2 * item2;
	}
}
double GENROSEvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 1.0;
	double item1, item2;
	for (i = 1; i < n; i++)
	{
		item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
		item2 = x[i] - 1.0;
		fx += item1*item1 + item2*item2;
	}
	return fx;
}
void GENROSEStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 1.0 / n + 1;
}
/*========= 33 ========== LIARWDH Function ============= 5000 ==============*/
double LIARWDHvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;
	double fx = (0.0);
	for (i = 0; i < n; i++)
	{
		item1 = 2.0*(x[i] * x[i] - x[0]);
		item2 = (x[i] - 1);
		fx += item1*item1 + item2*item2;
		g[i] = 8.0*item1*x[i] + 2.0*item2;
		g[0] -= 4.0*item1;

	}
	return fx;
}
void LIARWDHgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;
	for (i = 0; i < n; i++)
	{
		item1 = 2.0*(x[i] * x[i] - x[0]);
		item2 = (x[i] - 1);
		g[i] = 8.0*item1*x[i] + 2.0*item2;
		g[0] -= 4.0*item1;

	}
}
double LIARWDHvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;
	double fx = (0.0);
	for (i = 0; i < n; i++)
	{
		item1 = 2.0*(x[i] * x[i] - x[0]);
		item2 = (x[i] - 1);
		fx += item1*item1 + item2*item2;
	}
	return fx;
}
void LIARWDHStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 4.0;
}
/*========= 34 ========== MOREBV ============= 5000 ==============*/
double MOREBVvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double h = (1.0 / n);
	double element, item;
	double fx = (0.0);
	for (i = 0; i < n; i++)
		g[i] = 0;
	//first term
	element = x[1] + h + 1;
	item = 2 * x[1] - x[2] + (h*h*pow(element, 3)) / 2;
	fx += item*item;
	g[1] = 2 * item*(2.0 + h*h*1.5*pow(element, 2));
	g[2] = -2 * item;
	//last term
	element = x[n - 1] + (n - 1)*h + 1;
	item = 2 * x[n - 1] - x[n - 2] + (h*h*pow(element, 3)) / 2;
	fx += item*item;
	g[n - 1] = 4 * item*(2.0 + h*h*1.5*pow(element, 2));
	g[n - 2] = -2 * item;
	for (i = 2; i < n - 1; i++)
	{
		element = x[i] + i*h + 1;
		item = 2 * x[i] - x[i - 1] - x[i + 1] + (h*h*pow(element, 3)) / 2;
		fx += item*item;
		g[i - 1] -= 2 * item;
		g[i] += 2 * item*(2.0 + h*h*1.5*pow(element, 2));
		g[i + 1] -= 2 * item;
	}
	return fx;
}
void MOREBVgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double h = (1.0 / n);
	double element, item;
	for (i = 0; i < n; i++)
		g[i] = 0;
	//first term
	element = x[1] + h + 1;
	item = 2 * x[1] - x[2] + (h*h*pow(element, 3)) / 2;
	g[1] = 2 * item*(2.0 + h*h*1.5*pow(element, 2));
	g[2] = -2 * item;
	//last term
	element = x[n - 1] + (n - 1)*h + 1;
	item = 2 * x[n - 1] - x[n - 2] + (h*h*pow(element, 3)) / 2;
	g[n - 1] = 4 * item*(2.0 + h*h*1.5*pow(element, 2));
	g[n - 2] = -2 * item;
	for (i = 2; i < n - 1; i++)
	{
		element = x[i] + i*h + 1;
		item = 2 * x[i] - x[i - 1] - x[i + 1] + (h*h*pow(element, 3)) / 2;
		g[i - 1] -= 2 * item;
		g[i] += 2 * item*(2.0 + h*h*1.5*pow(element, 2));
		g[i + 1] -= 2 * item;
	}
}
double MOREBVvalue
(
	double *x,
	INT n
)
{
	INT i;
	double h = (1.0 / n);
	double element, item;
	double fx = (0.0);
	//first term
	element = x[1] + h + 1;
	item = 2 * x[1] - x[2] + (h*h*pow(element, 3)) / 2;
	fx += item*item;
	//last term
	element = x[n - 1] + (n - 1)*h + 1;
	item = 2 * x[n - 1] - x[n - 2] + (h*h*pow(element, 3)) / 2;
	fx += item*item;
	for (i = 2; i < n - 1; i++)
	{
		element = x[i] + i*h + 1;
		item = 2 * x[i] - x[i - 1] - x[i + 1] + (h*h*pow(element, 3)) / 2;
		fx += item*item;
	}
	return fx;
}
void MOREBVStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	double h = (1.0 / n);
	x[0] = 0.;
	for (i = 1; i< n; i++)
		x[i] = i*h*(i*h - 1.0);
}
/*===== 35 ======================= NONCVXU2 ============ 1000 ===============*/
double NONCVXU2valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	double fx = 0.0;
	for (i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2;
	for (j = 0; j < n; j++)
	{
		int i = (j + 1);
		item1 = (x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n]);
		item2 = x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n];
		fx += item1*item1 + 4 * cos(item2);
		g[j] += 2 * item1 - 4 * sin(item2);
		g[(3 * i - 2) % n] += 2 * item1 - 4 * sin(item2);
		g[(7 * i - 3) % n] += 2 * item1 - 4 * sin(item2);
	}
	return fx;
}
void NONCVXU2grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i, j;
	for (i = 0; i < n; i++)
		g[i] = 0;
	double item1, item2;
	for (j = 0; j < n; j++)
	{
		int i = (j + 1);
		item1 = (x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n]);
		item2 = x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n];
		g[j] += 2 * item1 - 4 * sin(item2);
		g[(3 * i - 2) % n] += 2 * item1 - 4 * sin(item2);
		g[(7 * i - 3) % n] += 2 * item1 - 4 * sin(item2);
	}
}
double NONCVXU2value
(
	double *x,
	INT n
)
{
	INT j;
	double fx = 0.0;
	double item1, item2;
	for (j = 0; j < n; j++)
	{
		int i = (j + 1);
		item1 = (x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n]);
		item2 = x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n];
		fx += item1*item1 + 4 * cos(item2);
	}
	return fx;
}
void NONCVXU2StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = i + 1;
}
/*========= 36 ========== NONDIA Function ============= 10000 ==============*/
double NONDIAvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1.0);
	double fx = (item*item);

	g[0] = 2.0*item;
	item = 10.0*(x[0] - x[0] * x[0]);
	fx += item*item;
	g[0] += (20.0 - 40.0*x[0])*item;

	for (i = 1; i < n; i++)
	{
		item = 10.0*(x[0] - x[i] * x[i]);
		fx += item*item;
		g[0] += 20.0*item;
		g[i] = -40.0*x[i] * item;
	}
	return fx;
}
void NONDIAgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1.0);

	g[0] = 2.0*item;
	item = 10.0*(x[0] - x[0] * x[0]);
	g[0] += (20.0 - 40.0*x[0])*item;

	for (i = 1; i < n; i++)
	{
		item = 10.0*(x[0] - x[i] * x[i]);
		g[0] += 20.0*item;
		g[i] = -40.0*x[i] * item;
	}
}double NONDIAvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1.0);
	double fx = (item*item);

	item = 10.0*(x[0] - x[0] * x[0]);
	fx += item*item;
	for (i = 1; i < n; i++)
	{
		item = 10.0*(x[0] - x[i] * x[i]);
		fx += item*item;
	}
	return fx;
}
void NONDIAStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i++)
		x[i] = -1.0;
}
/*========= 37 ========== NONDQUAR Function ============= 10000 ==============*/
double NONDQUARvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = (0.0);

	for (i = 0; i < n; i++)
		g[i] = 0.0;

	item = x[0] - x[1];
	fx += item*item;
	g[0] = 2.0 * item;
	g[1] = -2.0 * item;

	item = x[n - 2] + x[n - 1];
	fx += item*item;
	g[n - 2] = 2.0 * item;
	g[n - 1] = 2.0 * item;

	for (i = 0; i < n - 2; i++)
	{
		item = x[i] + x[i + 1] + x[n - 1];
		fx += pow(item, 4.0);
		g[i] += 4.0*pow(item, 3.0);
		g[i + 1] += 4.0*pow(item, 3.0);
		g[n - 1] += 4.0*pow(item, 3.0);
	}
	return fx;
}
void NONDQUARgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;

	for (i = 0; i < n; i++)
		g[i] = 0.0;

	item = x[0] - x[1];
	g[0] = 2.0 * item;
	g[1] = -2.0 * item;

	item = x[n - 2] + x[n - 1];
	g[n - 2] = 2.0 * item;
	g[n - 1] = 2.0 * item;

	for (i = 0; i < n - 2; i++)
	{
		item = x[i] + x[i + 1] + x[n - 1];
		g[i] += 4.0*pow(item, 3.0);
		g[i + 1] += 4.0*pow(item, 3.0);
		g[n - 1] += 4.0*pow(item, 3.0);
	}
}
double NONDQUARvalue
(
	double *x,
	INT n
)
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
void NONDQUARStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i += 2)
	{
		x[i] = 1.0;
		x[i + 1] = -1.0;
	}
}
/*========= 38 ========== PENALTY1 ============= 1000 ==============*/
double PENALTY1valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double tail = (0.0);
	double a = (1E-5);
	double fx = 0.0;

	for (i = 0; i < n; i++) {
		item = x[i] - 1;
		fx += a * item*item;
		g[i] = 2.0*a * item;
		tail += pow(x[i], 2);
	}
	fx += pow(tail - 0.25, 2);
	for (i = 0; i < n; i++)
		g[i] += 4 * (tail - 0.25)*x[i];
	return fx;
}
void PENALTY1grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double tail = (0.0);
	double a = (1E-5);
	double fx = 0.0;

	for (i = 0; i < n; i++) {
		item = x[i] - 1;
		g[i] = 2.0*a * item;
		tail += pow(x[i], 2);
	}
	for (i = 0; i < n; i++)
		g[i] += 4 * (tail - 0.25)*x[i];
}
double PENALTY1value
(
	double *x,
	INT n
)
{
	INT i;
	double item;
	double tail = (0.0);
	double a = (1E-5);
	double fx = 0.0;

	for (i = 0; i < n; i++) {
		item = x[i] - 1;
		fx += a * item*item;
		tail += pow(x[i], 2);
	}
	fx += pow(tail - 0.25, 2);
	return fx;
}

void PENALTY1StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i++)
		x[i] = i + 1;
}
/*========= 39 ========== PENALTY2 ============= 100 ==============*/
double PENALTY2valgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;
	double tail = (0.0);
	double a = (1E-5);
	for (i = 0; i < n; i++)
		g[i] = 0.0;
	item1 = x[0] - 0.2;
	double fx = (pow(item1, 2));
	g[0] = 2 * item1;
	tail += (n)*pow(x[0], 2);
	for (i = 1; i < n; i++)
	{
		item1 = exp(x[i] / 10) + exp(x[i - 1] / 10) - exp((i + 1) / 10.0) - exp(i / 10.0);
		item2 = exp(x[i] / 10) - exp(-1 / 10.0);
		tail += (n - i)*pow(x[i], 2);
		fx += a * (item1*item1 + item2*item2);
		g[i - 1] += 0.2* a * item1*exp(x[i - 1] / 10);
		g[i] += 0.2* a * exp(x[i] / 10)*(item1 + item2);
	}
	fx += (tail - 1)*(tail - 1);
	for (i = 0; i < n; i++)
		g[i] += 4 * (tail - 1)*(n - i)*x[i];
	return fx;
}
void PENALTY2grad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;
	double tail = (0.0);
	double a = (1E-5);
	for (i = 0; i < n; i++)
		g[i] = 0.0;
	item1 = x[0] - 0.2;
	g[0] = 2 * item1;
	tail += (n)*pow(x[0], 2);
	for (i = 1; i < n; i++)
	{
		item1 = exp(x[i] / 10) + exp(x[i - 1] / 10) - exp((i + 1) / 10.0) - exp(i / 10.0);
		item2 = exp(x[i] / 10) - exp(-1 / 10.0);
		tail += (n - i)*pow(x[i], 2);
		g[i - 1] += 0.2* a * item1*exp(x[i - 1] / 10);
		g[i] += 0.2* a * exp(x[i] / 10)*(item1 + item2);
	}
	for (i = 0; i < n; i++)
		g[i] += 4 * (tail - 1)*(n - i)*x[i];
}
double PENALTY2value
(
	double *x,
	INT n
)
{
	INT i;
	double item1, item2;
	double tail = (0.0);
	double a = (1E-5);
	item1 = x[0] - 0.2;
	double fx = (pow(item1, 2));
	tail += (n)*pow(x[0], 2);
	for (i = 1; i < n; i++)
	{
		item1 = exp(x[i] / 10) + exp(x[i - 1] / 10) - exp((i + 1) / 10.0) - exp(i / 10.0);
		item2 = exp(x[i] / 10) - exp(-1 / 10.0);
		tail += (n - i)*pow(x[i], 2);
		fx += a * (item1*item1 + item2*item2);
	}
	fx += (tail - 1)*(tail - 1);
	return fx;
}
void PENALTY2StartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i++)
		x[i] = 0.5;

}
/*========= 40 ========== POWER Function ======== 1000 ==========*/
double POWERvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = (0.0);

	for (i = 0; i < n; i++)
	{
		item = (i + 1)*x[i];
		fx += item*item;
		g[i] = 2.0*item*(i + 1);
	}
	return fx;
}
void POWERgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	for (i = 0; i < n; i++)
	{
		item = (i + 1)*x[i];
		g[i] = 2.0*item*(i + 1);
	}
}
double POWERvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = (0.0);

	for (i = 0; i < n; i++)
	{
		item = (i + 1)*x[i];
		fx += item*item;
	}
	return fx;
}
void POWERStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i++)
		x[i] = 1.0;

}
/*========= 41 ====== QARTC Function ======== 10000 =================*/
double QARTCvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = 0.0;

	for (i = 0; i < n; i++)
	{
		item = (x[i] - i - 1);
		fx += pow(item, 4.0);
		g[i] = 4 * item* item* item;
	}
	return fx;
}
void QARTCgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item;
	for (i = 0; i < n; i++)
	{
		item = (x[i] - i - 1);
		g[i] = 4 * item* item* item;
	}
}
double QARTCvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item;
	double fx = 0.0;

	for (i = 0; i < n; i++)
	{
		item = (x[i] - i - 1);
		fx += pow(item, 4.0);
	}
	return fx;
}
void QARTCStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i++)
		x[i] = 2.0;
}
/*===== 42 =============== SROSENBR ============= 10000 ================*/
double SROSENBRvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	for (i = 0; i < n; i += 2) {
		double t1 = 1.0 - x[i];
		double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
		g[i + 1] = 20.0 * t2;
		g[i] = -2.0 * (x[i] * g[i + 1] + t1);
		fx += t1 * t1 + t2 * t2;
	}
	return fx;
}
void SROSENBRgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i += 2) {
		double t1 = 1.0 - x[i];
		double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
		g[i + 1] = 20.0 * t2;
		g[i] = -2.0 * (x[i] * g[i + 1] + t1);
	}
}
double SROSENBRvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = 0.0;
	for (i = 0; i < n; i += 2) {
		double t1 = 1.0 - x[i];
		double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
		fx += t1 * t1 + t2 * t2;
	}
	return fx;
}
void SROSENBRStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i < n; i += 2)
	{
		x[i] = -1.2;
		x[i + 1] = 1.0;
	}
}
/*========= 43 ========== TRIDIA Function ============= 10000 ==============*/
double TRIDIAvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1.0);
	double fx = (item*item);
	for (i = 0; i < n; i++)
		g[i] = 0;

	g[0] = 2.0*item;
	for (i = 1; i < n; i++)
	{
		item = (2.0*x[i] - x[i - 1]);
		fx += (i + 1)* item*item;
		g[i] += 4.0*item*(i + 1);
		g[i - 1] -= 2.0*(i + 1)*item;
	}
	return fx;
}
void TRIDIAgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1.0);
	for (i = 0; i < n; i++)
		g[i] = 0;

	g[0] = 2.0*item;
	for (i = 1; i < n; i++)
	{
		item = (2.0*x[i] - x[i - 1]);
		g[i] += 4.0*item*(i + 1);
		g[i - 1] -= 2.0*(i + 1)*item;
	}
}
double TRIDIAvalue
(
	double *x,
	INT n
)
{
	INT i;
	double item = (x[0] - 1.0);
	double fx = (item*item);
	for (i = 1; i < n; i++)
	{
		item = (2.0*x[i] - x[i - 1]);
		fx += (i + 1)* item*item;
	}
	return fx;
}
void TRIDIAStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i++)
		x[i] = 1.0;
}
/*========= 44 =========  Woods ======== 10000 =====================*/
double WOODSvalgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double fx = (0.0);
	double item1, item2, item3, item4, item5, item6;

	for (i = 0; i < n; i += 4)
	{
		item1 = (x[i + 1] - pow(x[i], 2.0));
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - pow(x[i + 2], 2.0));
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		fx += (100 * pow(item1, 2.0) + pow(item2, 2.0) + 90 * pow(item3, 2.0)
			+ pow(item4, 2.0) + 10.0*(pow(item5, 2.0) + 0.1*pow(item6, 2.0)));
		g[i] = -400 * item1*x[i] - 2 * item2;
		g[i + 1] = 200 * item1 + 20.0*item5 + 0.2*item6;
		g[i + 2] = -360 * item3*x[i + 2] - 2.0*item4;
		g[i + 3] = 180 * item3 + 20.0*item5 - 0.2*item6;
	}
	return fx;
}
void WOODSgrad
(
	double *g,
	double *x,
	INT n
)
{
	INT i;
	double item1, item2, item3, item4, item5, item6;

	for (i = 0; i < n; i += 4)
	{
		item1 = (x[i + 1] - pow(x[i], 2.0));
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - pow(x[i + 2], 2.0));
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		g[i] = -400 * item1*x[i] - 2 * item2;
		g[i + 1] = 200 * item1 + 20.0*item5 + 0.2*item6;
		g[i + 2] = -360 * item3*x[i + 2] - 2.0*item4;
		g[i + 3] = 180 * item3 + 20.0*item5 - 0.2*item6;
	}
}
double WOODSvalue
(
	double *x,
	INT n
)
{
	INT i;
	double fx = (0.0);
	double item1, item2, item3, item4, item5, item6;

	for (i = 0; i < n; i += 4)
	{
		item1 = (x[i + 1] - pow(x[i], 2.0));
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - pow(x[i + 2], 2.0));
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		fx += (100 * pow(item1, 2.0) + pow(item2, 2.0) + 90 * pow(item3, 2.0)
			+ pow(item4, 2.0) + 10.0*(pow(item5, 2.0) + 0.1*pow(item6, 2.0)));
	}
	return fx;
}
void WOODSStartingGess
(
	double *x,
	INT n
)
{
	INT i;
	for (i = 0; i<n; i += 2)
	{
		x[i] = -3;
		x[i + 1] = -1;
	}
}
/*------------------  end of tests ------------------------*/
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
const INT TestNumber = 45;
struct testFunction test[45];

void makeTestCollection(struct testFunction* t)
{
	//1-ARWHEAD
	test[0].size = 5000;
	strcpy(test[0].name, "ARWHEAD");
	test[0].value = ARWHEADvalue;
	test[0].grad = ARWHEADgrad;
	test[0].valgrad = ARWHEADvalgrad;
	test[0].StartingGess = ARWHEADStartingGess;
	//2-BDQRTIC
	test[1].size = 5000;
	strcpy(test[1].name, "BDQRTIC");
	test[1].value = BDQRTICvalue;
	test[1].grad = BDQRTICgrad;
	test[1].valgrad = BDQRTICvalgrad;
	test[1].StartingGess = BDQRTICStartingGess;
	//3-BROYDN7D
	test[2].size = 5000;
	strcpy(test[2].name, "BROYDN7D");
	test[2].value = BROYDN7Dvalue;
	test[2].grad = BROYDN7Dgrad;
	test[2].valgrad = BROYDN7Dvalgrad;
	test[2].StartingGess = BROYDN7DStartingGess;
	//4-BRYBND 
	test[3].size = 5000;
	strcpy(test[3].name, "BRYBND");
	test[3].value = BRYBNDvalue;
	test[3].grad = BRYBNDgrad;
	test[3].valgrad = BRYBNDvalgrad;
	test[3].StartingGess = BRYBNDStartingGess;
	//5-CHAINWOO
	test[4].size = 4000;
	strcpy(test[4].name, "CHAINWOO");
	test[4].value = CHAINWOOvalue;
	test[4].grad = CHAINWOOgrad;
	test[4].valgrad = CHAINWOOvalgrad;
	test[4].StartingGess = CHAINWOOStartingGess;
	//6-COSINE
	test[5].size = 10000;
	strcpy(test[5].name, "COSINE");
	test[5].value = COSINEvalue;
	test[5].grad = COSINEgrad;
	test[5].valgrad = COSINEvalgrad;
	test[5].StartingGess = COSINEStartingGess;
	//7-CRAGGLVY
	test[6].size = 5000;
	strcpy(test[6].name, "CRAGGLVY");
	test[6].value = CRAGGLVYvalue;
	test[6].grad = CRAGGLVYgrad;
	test[6].valgrad = CRAGGLVYvalgrad;
	test[6].StartingGess = CRAGGLVYStartingGess;
	//8-CURLY10
	test[7].size = 500;
	strcpy(test[7].name, "CURLY10");
	test[7].value = CURLY10value;
	test[7].grad = CURLY10grad;
	test[7].valgrad = CURLY10valgrad;
	test[7].StartingGess = CURLY10StartingGess;
	//9-CURLY20
	test[8].size = 500;
	strcpy(test[8].name, "CURLY20");
	test[8].value = CURLY20value;
	test[8].grad = CURLY20grad;
	test[8].valgrad = CURLY20valgrad;
	test[8].StartingGess = CURLY20StartingGess;
	//10-CURLY30
	test[9].size = 500;
	strcpy(test[9].name, "CURLY30");
	test[9].value = CURLY30value;
	test[9].grad = CURLY30grad;
	test[9].valgrad = CURLY30valgrad;
	test[9].StartingGess = CURLY30StartingGess;
	//11-DIXMAANA
	test[10].size = 5000;
	strcpy(test[10].name, "DIXMAANA");
	test[10].value = DIXMAANAvalue;
	test[10].grad = DIXMAANAgrad;
	test[10].valgrad = DIXMAANAvalgrad;
	test[10].StartingGess = DIXMAANStartingGess;
	//12-DIXMAANB
	test[11].size = 5000;
	strcpy(test[11].name, "DIXMAANB");
	test[11].value = DIXMAANBvalue;
	test[11].grad = DIXMAANBgrad;
	test[11].valgrad = DIXMAANBvalgrad;
	test[11].StartingGess = DIXMAANStartingGess;
	//13-DIXMAANC
	test[12].size = 5000;
	strcpy(test[12].name, "DIXMAANC");
	test[12].value = DIXMAANCvalue;
	test[12].grad = DIXMAANCgrad;
	test[12].valgrad = DIXMAANCvalgrad;
	test[12].StartingGess = DIXMAANStartingGess;
	//14-DIXMAAND
	test[13].size = 5000;
	strcpy(test[13].name, "DIXMAAND");
	test[13].value = DIXMAANDvalue;
	test[13].grad = DIXMAANDgrad;
	test[13].valgrad = DIXMAANDvalgrad;
	test[13].StartingGess = DIXMAANStartingGess;
	//DIXMAANE
	test[14].size = 3000;
	strcpy(test[14].name, "DIXMAANE");
	test[14].value = DIXMAANEvalue;
	test[14].grad = DIXMAANEgrad;
	test[14].valgrad = DIXMAANEvalgrad;
	test[14].StartingGess = DIXMAANStartingGess;
	//DIXMAANF
	test[15].size = 3000;
	strcpy(test[15].name, "DIXMAANF");
	test[15].value = DIXMAANFvalue;
	test[15].grad = DIXMAANFgrad;
	test[15].valgrad = DIXMAANFvalgrad;
	test[15].StartingGess = DIXMAANStartingGess;
	//DIXMAANG
	test[16].size = 3000;
	strcpy(test[16].name, "DIXMAANG");
	test[16].value = DIXMAANGvalue;
	test[16].grad = DIXMAANGgrad;
	test[16].valgrad = DIXMAANGvalgrad;
	test[16].StartingGess = DIXMAANStartingGess;
	//DIXMAANH
	test[17].size = 3000;
	strcpy(test[17].name, "DIXMAANH");
	test[17].value = DIXMAANHvalue;
	test[17].grad = DIXMAANHgrad;
	test[17].valgrad = DIXMAANHvalgrad;
	test[17].StartingGess = DIXMAANStartingGess;
	//DIXMAANI
	test[18].size = 3000;
	strcpy(test[18].name, "DIXMAANI");
	test[18].value = DIXMAANIvalue;
	test[18].grad = DIXMAANIgrad;
	test[18].valgrad = DIXMAANIvalgrad;
	test[18].StartingGess = DIXMAANStartingGess;
	//DIXMAANJ
	test[19].size = 3000;
	strcpy(test[19].name, "DIXMAANJ");
	test[19].value = DIXMAANJvalue;
	test[19].grad = DIXMAANJgrad;
	test[19].valgrad = DIXMAANJvalgrad;
	test[19].StartingGess = DIXMAANStartingGess;
	//DIXMAANK
	test[20].size = 3000;
	strcpy(test[20].name, "DIXMAANK");
	test[20].value = DIXMAANKvalue;
	test[20].grad = DIXMAANKgrad;
	test[20].valgrad = DIXMAANKvalgrad;
	test[20].StartingGess = DIXMAANStartingGess;
	//DIXMAANL
	test[21].size = 3000;
	strcpy(test[21].name, "DIXMAANL");
	test[21].value = DIXMAANLvalue;
	test[21].grad = DIXMAANLgrad;
	test[21].valgrad = DIXMAANLvalgrad;
	test[21].StartingGess = DIXMAANStartingGess;
	//23-DIXON3DQ
	test[22].size = 5000;
	strcpy(test[22].name, "DIXON3DQ");
	test[22].value = DIXON3DQvalue;
	test[22].grad = DIXON3DQgrad;
	test[22].valgrad = DIXON3DQvalgrad;
	test[22].StartingGess = DIXON3DQStartingGess;
	//24-DQDRTIC
	test[23].size = 10000;
	strcpy(test[23].name, "DQDRTIC");
	test[23].value = DQDRTICvalue;
	test[23].grad = DQDRTICgrad;
	test[23].valgrad = DQDRTICvalgrad;
	test[23].StartingGess = DQDRTICStartingGess;
	//25-DQRTIC
	test[24].size =  100000;
	strcpy(test[24].name, "DQRTIC");
	test[24].value = DQRTICvalue;
	test[24].grad = DQRTICgrad;
	test[24].valgrad = DQRTICvalgrad;
	test[24].StartingGess = DQRTICStartingGess;
	//26-EDENSCH
	test[25].size = 2000;
	strcpy(test[25].name, "EDENSCH");
	test[25].value = EDENSCHvalue;
	test[25].grad = EDENSCHgrad;
	test[25].valgrad = EDENSCHvalgrad;
	test[25].StartingGess = EDENSCHStartingGess;
	//27-EG2
	test[26].size = 5000;
	strcpy(test[26].name, "EG2");
	test[26].value = EG2value;
	test[26].grad = EG2grad;
	test[26].valgrad = EG2valgrad;
	test[26].StartingGess = EG2StartingGess;
	//28-ENGVAL1
	test[27].size = 5000;
	strcpy(test[27].name, "ENGVAL1");
	test[27].value = ENGVAL1value;
	test[27].grad = ENGVAL1grad;
	test[27].valgrad = ENGVAL1valgrad;
	test[27].StartingGess = ENGVAL1StartingGess;
	//EXTROSNB
	test[28].size = 1000;
	strcpy(test[28].name, "EXTROSNB");
	test[28].value = EXTROSNBvalue;
	test[28].grad = EXTROSNBgrad;
	test[28].valgrad = EXTROSNBvalgrad;
	test[28].StartingGess = EXTROSNBStartingGess;
	//FLETCHR
	test[29].size = 1000;
	strcpy(test[29].name, "FLETCHR");
	test[29].value = FLETCHRvalue;
	test[29].grad = FLETCHRgrad;
	test[29].valgrad = FLETCHRvalgrad;
	test[29].StartingGess = FLETCHRStartingGess;
	//31-FREUROTH
	test[30].size = 5000;
	strcpy(test[30].name, "FREUROTH");
	test[30].value = FREUROTHvalue;
	test[30].grad = FREUROTHgrad;
	test[30].valgrad = FREUROTHvalgrad;
	test[30].StartingGess = FREUROTHStartingGess;
	//32-GENHUMPS
	test[31].size = 5000;
	strcpy(test[31].name, "GENHUMPS");
	test[31].value = GENHUMPSvalue;
	test[31].grad = GENHUMPSgrad;
	test[31].valgrad = GENHUMPSvalgrad;
	test[31].StartingGess = GENHUMPSStartingGess;
	//33-GENROSE
	test[32].size = 1000;
	strcpy(test[32].name, "GENROSE");
	test[32].value = GENROSEvalue;
	test[32].grad = GENROSEgrad;
	test[32].valgrad = GENROSEvalgrad;
	test[32].StartingGess = GENROSEStartingGess;
	//34-LIARWDH
	test[33].size = 5000;
	strcpy(test[33].name, "LIARWDH");
	test[33].value = LIARWDHvalue;
	test[33].grad = LIARWDHgrad;
	test[33].valgrad = LIARWDHvalgrad;
	test[33].StartingGess = LIARWDHStartingGess;
	//35-MOREBV
	test[34].size = 5000;
	strcpy(test[34].name, "MOREBV");
	test[34].value = MOREBVvalue;
	test[34].grad = MOREBVgrad;
	test[34].valgrad = MOREBVvalgrad;
	test[34].StartingGess = MOREBVStartingGess;
	//36-NONCVXU2
	test[35].size = 1000;
	strcpy(test[35].name, "NONCVXU2");
	test[35].value = NONCVXU2value;
	test[35].grad = NONCVXU2grad;
	test[35].valgrad = NONCVXU2valgrad;
	test[35].StartingGess = NONCVXU2StartingGess;
	//37-NONDIA
	test[36].size = 10000;
	strcpy(test[36].name, "NONDIA");
	test[36].value = NONDIAvalue;
	test[36].grad = NONDIAgrad;
	test[36].valgrad = NONDIAvalgrad;
	test[36].StartingGess = NONDIAStartingGess;
	//38-NONDQUAR
	test[37].size = 10000;
	strcpy(test[37].name, "NONDQUAR");
	test[37].value = NONDQUARvalue;
	test[37].grad = NONDQUARgrad;
	test[37].valgrad = NONDQUARvalgrad;
	test[37].StartingGess = NONDQUARStartingGess;
	//39-PENALTY1
	test[38].size = 1000;
	strcpy(test[38].name, "PENALTY1");
	test[38].value = PENALTY1value;
	test[38].grad = PENALTY1grad;
	test[38].valgrad = PENALTY1valgrad;
	test[38].StartingGess = PENALTY1StartingGess;
	//40-PENALTY2
	test[39].size = 100;
	strcpy(test[39].name, "PENALTY2");
	test[39].value = PENALTY2value;
	test[39].grad = PENALTY2grad;
	test[39].valgrad = PENALTY2valgrad;
	test[39].StartingGess = PENALTY2StartingGess;
	//41-POWER
	test[40].size = 1000;
	strcpy(test[40].name, "POWER");
	test[40].value = POWERvalue;
	test[40].grad = POWERgrad;
	test[40].valgrad = POWERvalgrad;
	test[40].StartingGess = POWERStartingGess;
	//42-QARTC
	test[41].size = 10000;
	strcpy(test[41].name, "QARTC");
	test[41].value = QARTCvalue;
	test[41].grad = QARTCgrad;
	test[41].valgrad = QARTCvalgrad;
	test[41].StartingGess = QARTCStartingGess;
	//43-SROSENBR
	test[42].size = 10000;
	strcpy(test[42].name, "SROSENBR");
	test[42].value = SROSENBRvalue;
	test[42].grad = SROSENBRgrad;
	test[42].valgrad = SROSENBRvalgrad;
	test[42].StartingGess = SROSENBRStartingGess;
	//44-TRIDIA
	test[43].size = 10000;
	strcpy(test[43].name, "TRIDIA");
	test[43].value = TRIDIAvalue;
	test[43].grad = TRIDIAgrad;
	test[43].valgrad = TRIDIAvalgrad;
	test[43].StartingGess = TRIDIAStartingGess;
	//45-WOODS
	test[44].size = 10000;
	strcpy(test[44].name, "WOODS");
	test[44].value = WOODSvalue;
	test[44].grad = WOODSgrad;
	test[44].valgrad = WOODSvalgrad;
	test[44].StartingGess = WOODSStartingGess;
}
void iSort(double* a, int size)
{
	if (size>1)
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
