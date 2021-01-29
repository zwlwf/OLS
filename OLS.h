#ifndef _OLS_H
#define _OLS_H
#include <vector>
#include "Param.h"

class OLS {
	std::vector<double> t;
	std::vector<double> y;
	Param m_beta; // store the final model
	public:
	// read training data from file
	void readData(const char* fname);
	void setData(double*t, double *y, int n);
	OLS() {}
	double E(const Param&, void *args ); // error function
	void dE(const Param&, double*v, void*args); // gradE
	void ddE(const Param&, double**v, void* args); // Hessian E
	Param sd(Param beta0, int maxItr, double tol, double alpha, void*args);
	Param newton(Param beta0, int maxItr, double tol, double alpha, void* args);
	Param trust_region(Param beta0, int maxItr, double tol, void*args);
};

#endif
