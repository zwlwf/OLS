#include "nonLinearEq.h"
#include <bits/stdc++.h>
#include "OLS.h"

using namespace std;

int n;
const int maxn = 1e6+5;
double t[maxn];
double y[maxn];
double theta[maxn];

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

void f(double *x, double *y, void* p) {
	OLS* ls = (OLS*) p;
	double omega = 9.8;
	ls->dE( Param(x[0], x[1], x[2]), y, (void*)&omega );
}

void df(double *x, double **H, void*p ) {
	OLS* ls = (OLS*) p;
	double omega = 9.8;
	ls->ddE( Param(x[0], x[1], x[2]), H, (void*)&omega );
}

int main() {
	readData("001.sec");
	OLS ls;
	ls.setData(t, y, n);
	nonLinearEq eq(3, f, df, (void*)&ls);
	double x0[3] = {0.005*9, 0.1, 0.1} ;
	double x[3];
	eq.sd( x0, 100, 1.e-6, 1.0, x);
	cout<<"x = "<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
	//eq.newton( &x0, 100, 1.e-6, 1.0, &x, 1);
	//cout<<"x = "<<x<<endl;
	return 0;
}

