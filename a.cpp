#include <bits/stdc++.h>
using namespace std;
#include "OLS.h"

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

int main() {
	readData("001.sec");
	OLS ls;
	ls.setData(t, y, n);
	double omega = 9.821;
	Param beta0(0.005*9, 0.1, 0.1);
	//ls.sd(beta0, 100, 1.0e-5, 1.0e-3, (void*) &omega);
	//ls.newton(beta0, 100, 1.0e-5, 1, (void*) &omega);
	ls.trust_region(beta0, 100, 1.0e-5, (void*) &omega);
	return 0;
}
