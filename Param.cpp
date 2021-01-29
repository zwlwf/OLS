#include "Param.h"
#include <bits/stdc++.h>
using namespace std;
const double pi = 3.141692653;

//---------- common APIs ------------------
void Param::operator=(const Param& o) {
	memcpy(this, &o, sizeof(o));
}

Param Param::operator*(const double alpha) {
	Param ans;
	for(int i=0; i<m; i++) ans.vec()[i] = this->vec()[i]*alpha;
	return ans;
}

void Param::operator+=(const Param& o) {
	for(int i=0; i<m; i++) this->vec()[i] += o.vec()[i];
}

Param Param::operator+(const Param& o) {
	Param ans;
	for(int i=0; i<m; i++) ans.vec()[i] = this->vec()[i] + o.vec()[i];
	return ans;
}

Param Param::operator-(const Param& o) {
	Param ans;
	for(int i=0; i<m; i++) ans.vec()[i] = this->vec()[i] - o.vec()[i];
	return ans;
}

//---------- Special APIs ------------------
double Param::h( double t, void* args) const{
	double omega = *(double*) args;
	return exp(-lambda*t) *( A*cos(omega*t) + B*sin(omega*t));
}

void Param::dh(double t, double*v, void* args) const {
	double omega = *(double*) args;
	v[0] = -t*exp(-lambda*t) *( A*cos(omega*t) + B*sin(omega*t));
	v[1] = exp(-lambda*t)*cos(omega*t);
	v[2] = exp(-lambda*t)*sin(omega*t);
}

void Param::ddh(double t, double**v, void* args) const{
	double omega = *(double*) args;
	double yhat = this->h(t, args);
	v[0][0] = t*t*exp(-lambda*t)*( A*cos(omega*t) + B*sin(omega*t) );
	v[0][1] = v[1][0] = -t*exp(-lambda*t)*cos(omega*t);
	v[0][2] = v[2][0] = -t*exp(-lambda*t)*sin(omega*t);
	v[1][1] = 0; 
	v[1][2] = v[2][1] = 0; 
	v[2][2] = 0;
}

void Param::disp() {
	double omega = 9.821;
	double xi = lambda/sqrt( lambda*lambda + omega*omega);
	double w = omega/sqrt(1-xi*xi);
	cout<<"freq = "<< w/2/pi << ", " <<"xi = "<< xi <<", " <<"A = "<<A<<", " <<"B = "<<B<<endl;
}
