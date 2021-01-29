#ifndef _NONLINEAREQ_H
#define _NONLINEAREQ_H

typedef void (*Fun)(double *x, double *y, void* param);
typedef void (*DFun)(double *x, double**H, void* param);
class nonLinearEq {
	int n;
	void *param; // latentavariable
	Fun f;
	DFun df;
	double norm(double *y);
	public:
	// newton's method
	void newton(double *x0, int maxItr, double tol, double alpha, //input
				   double *x);
	// steepest descent method
	void sd(double *x0, int maxItr, double tol, double alpha,//input
				   double *x);
	nonLinearEq(int _n, Fun _f, DFun _df, void *_param) {
		n = _n;
		param = _param;
		f = _f;
		df = _df;
	}
};
#endif
