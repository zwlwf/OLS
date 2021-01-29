#include "OLS.h"
#include <bits/stdc++.h>

using namespace std;

double OLS::E( const Param& beta, void* args) {
	double res = 0;
	for(int i=0; i<t.size(); i++) {
		res += pow( y[i] - beta.h(t[i], args), 2);
	}
	return 0.5*res;
}

void OLS::dE( const Param& beta, double *v, void *args) {
	int m = beta.m;
	int n = t.size();
	double *tv = new double[m];
	for(int i=0; i<m; i++) v[i] = 0.0;
	for(int i=0; i<n; i++) {
		double yhat = beta.h(t[i], args);
		double res = y[i] - yhat;
		beta.dh(t[i], tv, args);
		for(int j=0; j<m; j++) v[j] -= res*tv[j];
	}
	delete []tv;
}

void OLS::ddE(const Param& beta, double**A, void*args) {
	int m = beta.m;
	int n = t.size();
	for(int i=0; i<m; i++) 
		for(int j=0; j<m; j++) A[i][j] = 0.0;

	double *v = new double[m];
	double **tA = new double*[m];
	double yhat;
	for(int j=0; j<m; j++) tA[j] = new double[m];
	for(int i=0; i<n; i++) {
		yhat = beta.h(t[i], args);
		beta.dh(t[i], v, args);
		beta.ddh(t[i], tA, args);
		for(int j=0; j<m; j++) 
			for(int k=0; k<m; k++) {
				A[j][k] += v[j]*v[k] - (y[i] - yhat)*tA[j][k];
			}
	}
	for(int j=0; j<m; j++) delete []tA[j];
	delete []tA;
	delete []v;
}

void OLS::setData(double*t ,double* y, int n) {
	this->t = vector<double>(n);
	this->y = vector<double>(n);
	for(int i=0; i<n; i++) {
		this->t[i] = t[i];
		this->y[i] = y[i];
	}
}

#define ZERO 1e-10
static void gausselim(double **A, double *b, int n) {
	for(int i=0; i<n; i++) {
		int maxI = i;
		double maxV = fabs(A[i][i]);
		for(int j=i+1; j<n; j++) {
			if(fabs(A[j][i]) > maxV ) {
				maxV = fabs(A[j][i]);
				maxI = j;
			}
		}
		if( maxV <= ZERO )
		{
			cout<<"Matrix is singular"<<endl;
			return ;
		}
		if(maxI != i) { // swap i, maxI
			for(int k=i; k<n; k++) swap(A[i][k], A[maxI][k]);
			swap(b[i], b[maxI]);
		}

		for(int j=i+1; j<n; j++) {
			double lambda = A[j][i]/A[i][i];
			for(int k=i; k<n; k++) {
				A[j][k] = A[j][k] - A[i][k]*lambda;
			}
			b[j] = b[j] - b[i]*lambda;
		}
	}
	// back trace
	for(int i=n-1; i>=0; i--) {
		double tmp = 0.0;
		for(int j=i+1; j<n; j++) tmp+= A[i][j]*b[j];
		b[i] = (b[i]-tmp)/A[i][i];
	}
}

// alpha is the step size
Param OLS::sd(Param beta0, int maxItr, double tol, double alpha, void* args) {
	double err = 1.5*tol;
	int itr = 0;
	m_beta = beta0;

	double *v=  new double[m_beta.m];
	while(err> tol && itr <maxItr ) {
		itr++;
		err = 0.0;
		dE(m_beta, v, args);
		for(int i=0; i<m_beta.m; i++){
		   	err += v[i]*v[i];
			m_beta.vec()[i] -= alpha*v[i];
		}
		err = sqrt(err);
		cout<<"iter = "<<itr<<", err= "<<err<<endl;
	}
	cout<<"iter = "<<itr<<", err= "<<err<<endl;
	m_beta.disp();
	delete []v;
	return m_beta;
}

Param OLS::newton(Param beta0, int maxItr, double tol, double alpha, void* args) {
	double err = 2*tol;
	int itr = 0;
	m_beta = beta0;
	int m = m_beta.m;
	int n = t.size();
	double **A = new double*[m];
	for(int i=0; i<m; i++) A[i] = new double[m];
	double *rhs = new double[m];
	double* tp = m_beta.vec();

	while( err>tol && itr<maxItr) {
		itr++;
		err = 0.0;
		ddE(m_beta, A, args);
		dE(m_beta, rhs, args);
		gausselim(A, rhs, m);
		for(int i=0; i<m; i++) {
			err += rhs[i]*rhs[i];
			tp[i] -= alpha*rhs[i];
		}
		err = sqrt(err);
		cout<<"iter = "<<itr<<", err= "<<err<<endl;
	}
		
	for(int i=0; i<m; i++) delete[] A[i];
	delete []A;
	delete []rhs;
	m_beta.disp();
	return m_beta;
}

Param OLS::trust_region(Param beta0, int maxItr, double tol, void*args) {
	// parameters for Trust region
	double t1 = 0.25;
	double t2 = 1.5;
	double eta0 = 1.0e-6;
	double eta1 = 0.25;
	double eta2 = 0.75;
	double DeltaMax = 10;

	double err = 2*tol;
	double res =2.0*tol;
	int itr = 0;
	m_beta = beta0;
	int m = m_beta.m;
	int n = t.size();
	double **A = new double*[m];
	for(int i=0; i<m; i++) A[i] = new double[m];

	Param dp;
	double *rhs = dp.vec();
	double* tp = m_beta.vec();

	vector<double> grad(m);
	vector<vector<double>> B(m, vector<double>(m));
	const double *s;
	double rho;
	double Delta = 1;

	while( itr<maxItr) {
		// cout<<"Delta = "<<Delta<<endl;
		err = 0.0;
		ddE(m_beta, A, args);
		dE(m_beta, rhs, args);
		// bak data before solve
		for(int i=0; i<m; i++) {
			grad[i] = rhs[i];
			err += rhs[i]*rhs[i];
			for(int j=0; j<m; j++)
				B[i][j] = A[i][j];
		}
		err = sqrt(err);
		if(err<tol) break;
#ifdef DEBUG
		//cout<<"iter = "<<itr<<", err= "<<err<<endl;
		m_beta.disp();
#endif
		gausselim(A, rhs, m); // solve with change the dp
		double norm_r = 0.0;
		for(int i=0; i<m; i++) norm_r += rhs[i]*rhs[i];
		norm_r = sqrt(norm_r);
		if(norm_r > Delta) {
			for(int i=0; i<m; i++) {
				rhs[i] = rhs[i]/norm_r*Delta;
			}
		}

		do {
			double nume = E(m_beta-dp, args) - E(m_beta, args);
			double deno = 0;
			s = rhs;
			for(int i=0; i<m; i++) {
				deno-= grad[i]*s[i];
				for(int j=0; j<m; j++) {
					deno += 0.5*s[i]*B[i][j]*s[j];
				}
			}
			rho = nume/deno;
			//cout<<nume<<" / "<<deno<<" "<<rho<<endl;
			if(rho > eta1 || rho<eta0) {
				for(int i=0; i<m; i++) {
					tp[i] -= rhs[i];
				}
				itr++;
			}

			if(rho >= eta2) {
				Delta = min(t2*Delta, DeltaMax);
				err = 1.0;
			} else if( rho <= eta1) {
				Delta *= t1;
				err = 1.0;
			} 
		} while(0);
			
	}
		
	for(int i=0; i<m; i++) delete[] A[i];
	delete []A;
	m_beta.disp();
	return m_beta;
}

