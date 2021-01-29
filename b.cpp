#include <bits/stdc++.h>

using namespace std;

int n;
const int maxn = 1e6+5;
double t[maxn];
double y[maxn];
double theta[maxn];
double res[maxn];

const double pi = 3.141692653;
struct Param {
	double lambda;
	double A;
	double B;
	double omega;
	void disp() {
		double xi = lambda/sqrt( lambda*lambda + omega*omega);
		double w = omega/sqrt(1-xi*xi);
		cout<<"freq = "<< w/2/pi << endl
			<<"xi = "<< xi <<endl
			<<"A = "<<A<<endl
			<<"B = "<<B<<endl;
	}
	double predict( double ti ) {
		return exp(-lambda*ti) *( A*cos(omega*ti) + B*sin(omega*ti));
	}

	// $\partial yhat/ \partial beta
	void df(double ti, double v[4]) {
		v[0] = -ti*exp(-lambda*ti) *( A*cos(omega*ti) + B*sin(omega*ti));
		v[1] = exp(-lambda*ti)*cos(omega*ti);
		v[2] = exp(-lambda*ti)*sin(omega*ti);
		v[3] = exp(-lambda*ti)*( -A*ti*sin(omega*ti) + B*ti*cos(omega*ti) ) ;
	}

	void d2f2(double ti, double v[4][4]) {
		double yhat = predict(ti);
		v[0][0] = ti*ti*exp(-lambda*ti)*( A*cos(omega*ti) + B*sin(omega*ti) );
		v[0][1] = v[1][0] = -ti*exp(-lambda*ti)*cos(omega*ti);
		v[0][2] = v[2][0] = -ti*exp(-lambda*ti)*sin(omega*ti);
		v[0][3] = v[3][0] = -ti*ti*exp(-lambda*ti)*( -A*sin(omega*ti) + B*cos(omega*ti) );
		v[1][1] = 0; 
		v[1][2] = v[2][1] = 0; 
		v[1][3] = v[3][1] = -ti*exp(-lambda*ti)*sin(omega*ti);
		v[2][2] = 0;
		v[2][3] = v[3][2] = ti*exp(-lambda*ti)*cos(omega*ti);
		v[3][3] = -ti*ti*exp(-lambda*ti)*(A*cos(omega*ti) + B*sin(omega*ti));
	}
};

void readData( const char* fname) {
	ifstream is(fname);	
	if(!is) {
		cout<<"import data failed!"<<endl;
		exit(-1);
	}
	string l;
	while( getline(is, l) ) {
		if(l.size()> 4 && l.substr(0,4) == "Time") break;
	}
	int tot = 0;
	while( is>>t[tot]>>y[tot]>>theta[tot] ) {
		tot++;
	}
	n = tot;
	is.close();
}

// training model
Param sd(double*t, double *y,  Param beta0, int maxItr, double tol) {
	double err = 1.5*tol;
	int itr = 0;
	double alpha = 1.0e-3; 
	Param p = beta0;
	Param dp;

	double res_norm;
	double v[4];
	while(err> tol && itr <maxItr ) {
		itr++;
		err = 0.0;
		dp.lambda = dp.A = dp.B = dp.omega = 0.0;
		res_norm = 0.0;
		for(int i=0; i<n; i++) {
			double yhat = p.predict(t[i]);
			res[i] = y[i] - yhat;
			res_norm += res[i]*res[i];
			p.df(t[i], v);
			dp.lambda -= res[i]*v[0];
			dp.A -= res[i]*v[1];
			dp.B -= res[i]*v[2];
			dp.omega -= res[i]*v[3];
		}

		p.lambda -= alpha*dp.lambda;
		p.A -= alpha*dp.A;
		p.B -= alpha*dp.B;
		p.omega -=alpha*dp.omega;
		err = dp.lambda*dp.lambda + dp.A*dp.A + dp.B*dp.B + dp.omega*dp.omega;
		err = sqrt(err);
		cout<<"iter = "<<itr<<", err= "<<err<<", res_norm= "<<res_norm<<endl;
	}
	cout<<"iter = "<<itr<<", err= "<<err<<endl;
	cout<<"res_norm = "<<res_norm<<endl;
	p.disp();
	return p;
}

void gausselim(double **A, double *b, int n) {
	for(int i=0; i<n; i++) {
		int maxI = i;
		double maxV = fabs(A[i][i]);
		for(int j=i+1; j<n; j++) {
			if(fabs(A[j][i]) > maxV ) {
				maxV = fabs(A[j][i]);
				maxI = j;
			}
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

// solve A x = rhs, and save the result x in rhs with gauss
void invArhs(double A[4][4], Param& rhs) {
	double *b = (double*) &rhs;
	auto print = [&] () {
		for(int ii=0; ii<4; ii++) {
			for(int jj=0; jj<4; jj++)
				cout<<A[ii][jj]<<" ";
			cout<<b[ii]<<endl;
		}
		cout<<"-------------- "<<endl;
	};

	int n = 4;
	double *Ain[4];
	for(int i=0; i<4; i++) Ain[i] = A[i];
	gausselim(Ain, b, 4);
}

// training model with newton's method
Param newton(double*tt, double *yy,  Param beta0, int maxItr, double tol) {
	double err = 2.5*tol;
	int itr = 0;
	double alpha = 0.1; 
	Param p;
	memcpy(&p, &beta0, sizeof(beta0));
	Param dp;

	int now = 0;
	double res_norm;
	double A[4][4];
	double dh[4];
	double ddh[4][4];
	while(err> tol && itr <maxItr ) {
		itr++;
		err = 0.0;
		for(int i=0; i<4; i++) for(int j=0; j<4; j++) A[i][j] = 0.0;
		dp.lambda = dp.A = dp.B = dp.omega = 0.0;
		res_norm = 0.0;

		for(int i=0; i<n; i++) {
			double yhat = p.predict(tt[i]);
			res[i] = yy[i] - yhat;
			res_norm += res[i]*res[i];
			p.df(tt[i], dh);
			p.d2f2(tt[i], ddh);
			dp.lambda -= res[i]*dh[0];
			dp.A -= res[i]*dh[1];
			dp.B -= res[i]*dh[2];
			dp.omega -= res[i]*dh[3];
			for(int k=0; k<4; k++)
				for(int j=0; j<4; j++) {
					A[k][j] += dh[k]*dh[j] - res[i]*ddh[k][j];
				}
		}

		invArhs(A, dp);
		p.lambda -= alpha*dp.lambda;
		p.A -= alpha*dp.A;
		p.B -= alpha*dp.B;
		p.omega -=alpha*dp.omega;
		err = dp.lambda*dp.lambda + dp.A*dp.A + dp.B*dp.B + dp.omega*dp.omega;
		err = sqrt(err);
		cout<<"iter = "<<itr<<", err= "<<err<<", res_norm= "<<res_norm<<endl;
	}
	cout<<"iter = "<<itr<<", err= "<<err<<endl;
	cout<<"res_norm = "<<res_norm<<endl;
	p.disp();
	return p;
}

// training model with trust region method
Param trust_region(double*tt, double *yy,  Param beta0, int maxItr, double tol) {
	double err = 2.5*tol;
	int itr = 0;
	double alpha = 1;  // init step size
	Param p;
	memcpy(&p, &beta0, sizeof(beta0));
	Param dp;

	double res_norm;
	double A[4][4];
	double dh[4];
	double ddh[4][4];
	while(err> tol && itr <maxItr ) {
		itr++;
		err = 0.0;
		for(int i=0; i<4; i++) for(int j=0; j<4; j++) A[i][j] = 0.0;
		dp.lambda = dp.A = dp.B = dp.omega = 0.0;
		res_norm = 0.0;

		for(int i=0; i<n; i++) {
			double yhat = p.predict(tt[i]);
			res[i] = yy[i] - yhat;
			res_norm += res[i]*res[i];
			p.df(tt[i], dh);
			p.d2f2(tt[i], ddh);
			dp.lambda -= res[i]*dh[0];
			dp.A -= res[i]*dh[1];
			dp.B -= res[i]*dh[2];
			dp.omega -= res[i]*dh[3];
			for(int k=0; k<4; k++)
				for(int j=0; j<4; j++) {
					A[k][j] += dh[k]*dh[j] - res[i]*ddh[k][j];
				}
		}

		invArhs(A, dp);
		p.lambda -= alpha*dp.lambda;
		p.A -= alpha*dp.A;
		p.B -= alpha*dp.B;
		p.omega -=alpha*dp.omega;
		err = dp.lambda*dp.lambda + dp.A*dp.A + dp.B*dp.B + dp.omega*dp.omega;
		err = sqrt(err);
		cout<<"iter = "<<itr<<", err= "<<err<<", res_norm= "<<res_norm<<endl;
	}
	cout<<"iter = "<<itr<<", err= "<<err<<endl;
	cout<<"res_norm = "<<res_norm<<endl;
	p.disp();
	return p;
}

// solve with fixed omega
Param newton2(double*tt, double *yy,  Param beta0, int maxItr, double tol) {
	double err = 2.5*tol;
	int itr = 0;
	double alpha = 0.1;  // init step size
	Param p;
	memcpy(&p, &beta0, sizeof(beta0));
	Param dp;

	double res_norm;
	double A[4][4];
	double dh[4];
	double ddh[4][4];
	while(err> tol && itr <maxItr ) {
		itr++;
		err = 0.0;
		for(int i=0; i<4; i++) for(int j=0; j<4; j++) A[i][j] = 0.0;
		dp.lambda = dp.A = dp.B = dp.omega = 0.0;
		res_norm = 0.0;

		for(int i=0; i<n; i++) {
			double yhat = p.predict(tt[i]);
			res[i] = yy[i] - yhat;
			res_norm += res[i]*res[i];
			p.df(tt[i], dh);
			p.d2f2(tt[i], ddh);
			dp.lambda -= res[i]*dh[0];
			dp.A -= res[i]*dh[1];
			dp.B -= res[i]*dh[2];
			dp.omega -= res[i]*dh[3];
			for(int k=0; k<4; k++)
				for(int j=0; j<4; j++) {
					A[k][j] += dh[k]*dh[j] - res[i]*ddh[k][j];
				}
		}

		//invArhs(A, dp);
		double *Ain[3];
		for(int i=0; i<3; i++) Ain[i] = A[i];
		gausselim(Ain, (double*) &dp, 3);
		p.lambda -= alpha*dp.lambda;
		p.A -= alpha*dp.A;
		p.B -= alpha*dp.B;
		//p.omega -=alpha*dp.omega;
		err = dp.lambda*dp.lambda + dp.A*dp.A + dp.B*dp.B ;
		err = sqrt(err);
		cout<<"iter = "<<itr<<", err= "<<err<<", res_norm= "<<res_norm<<endl;
	}
	cout<<"iter = "<<itr<<", err= "<<err<<endl;
	cout<<"res_norm = "<<res_norm<<endl;
	p.disp();
	return p;
}
int main() {
	readData("001.sec");
	// curve fit with steep-descent method
	Param h_param0{0.005*9, 0.1, 0.1,9.821};
	Param theta_param0{0.005*9, 0.1, 0,20.27};

	//Param ph = sd(t, y, h_param0, 100, 1.0e-5 );
	//Param pt = sd(t, theta, theta_param0, 1000, 1.0e-5 );
	Param ph = newton(t, y, h_param0, 100, 1.0e-5 );
	//Param pt = newton(t, theta, theta_param0, 100, 1.0e-5 );
	ofstream os("out");
	for(int i=0; i<n; i++) {
		os<<t[i]<<" "<<ph.predict(t[i])<<endl;
	}
	os.close();
	
	return 0;
}
