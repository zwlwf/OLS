#include "nonLinearEq.h"
#include <bits/stdc++.h>

double nonLinearEq::norm(double *y) {
	double ans = 0.0;
	for(int i=0; i<n; i++) ans += y[i]*y[i];
	return std::sqrt(ans);
}

void nonLinearEq::sd(double *x0, int maxItr, double tol, double alpha, //input
				   double *x) { // output
	memcpy( x, x0, n*sizeof(double) );
	int itr = 0;
	double *res = new double[n];
	while( itr < maxItr ) {
		itr++;
		f(x, res, param);
		std::cout<<"itr = "<<itr<<",err = "<<norm(res)<<std::endl;
		if( norm(res) < tol ) break;
		for(int i=0; i<n; i++)
			x[i] += alpha*res[i]; 
	}
	delete []res;
}

static void gausselim(double **A_, double *b, double *x, int n) {
	double** A = new double*[n];
	// copy matrix A = [A_ | b], not to pollute the original data
	for(int i=0; i<n; i++) A[i] = new double[n+1];
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) A[i][j] = A_[i][j];
		A[i][n] = b[i];
	}
	for(int i=0; i<n-1; i++) {
		int maxI = i;
		double maxV = fabs(A[i][i]);
		for(int j=i+1; j<n; j++) 
			if( fabs(A[j][i]) > maxV ) {
				maxV = fabs(A[j][i]);
				maxI = j;
			}
		if(maxI != i) {
			for(int j=i; j<n+1; j++)
				std::swap(A[i][j], A[maxI][j]);
		}
		for(int j=i+1; j<n; j++) {
			double la = A[j][i]/A[i][i];
			for(int k=i; k<n+1; k++)
				A[j][k] = A[j][k] - A[i][k]*la;
		}
	}
	// back substitution
	for(int i=n-1; i>=0; i--) {
		for(int j=i+1; j<n; j++) A[i][n] -= A[i][j]*x[j];
		x[i] = A[i][n]/A[i][i];
	}
	for(int i=0; i<n; i++) delete [] A[i];
	delete []A;
}

void nonLinearEq::newton(double *x0, int maxItr, double tol, double alpha, //input
				   double *x) {
	memcpy( x, x0, n*sizeof(double) );
	int itr = 0;
	double *res = new double[n];
	double *dx = new double[n];
	double **A = new double*[n];
	for(int i=0; i<n; i++) A[i] = new double[n];
	while( itr++ < maxItr ) {
		f(x, res, param);
		df(x, A, param);
		gausselim(A, res, dx, 1);
		if ( norm(res)<tol && norm(dx)<tol ) break;
		for(int i=0; i<n; i++)
			x[i] -= dx[i]*alpha;;
	}
	delete []dx;
	delete []res;
	for(int i=0; i<n; i++) delete []A[i];
	delete []A;
}

