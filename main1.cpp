#include "nonLinearEq.h"
#include <bits/stdc++.h>
#include "OLS.h"

using namespace std;

struct Param {
	double M;
	double e;
};

void f(double*x, double *y, void* param)
{
	Param p = *(Param*) param;
	y[0] = p.M + p.e*sin(x[0]) - x[0];
}

void df(double*x, double **H, void* param)
{
	Param p = *(Param*) param;
	H[0][0] = p.e*cos(x[0]) - 1;
}

int main() {
	Param p = {1, 0.5};
	nonLinearEq eq(1, f, df, (void*)&p);
	double x0 = 1;
	double x;
	eq.sd( &x0, 100, 1.e-6, 1.0, &x, 1);
	cout<<"x = "<<x<<endl;
	eq.newton( &x0, 100, 1.e-6, 1.0, &x, 1);
	cout<<"x = "<<x<<endl;
	return 0;
}

