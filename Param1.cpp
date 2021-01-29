#include "Param.h"
#include <bits/stdc++.h>
using namespace std;
const double pi = 3.141692653;

void Param::operator=(const Param& o) {
	lambda = o.lambda;
	A = o.A;
	B = o.B;
	omega = o.omega;
}

Param Param::operator*(const double alpha) {
	return Param( lambda*alpha, A*alpha, B*alpha, omega*alpha );
}

void Param::operator+=(const Param& o) {
	lambda += o.lambda;
	A += o.A;
	B += o.B;
	omega += o.omega;
}

Param Param::operator+(const Param& o) {
	Param ans;
	ans.lambda = lambda + o.lambda;
	ans.A = A + o.A;
	ans.B = B + o.B;
	ans.omega = omega + o.omega;
	return ans;
}

Param Param::operator-(const Param& o) {
	Param ans;
	ans.lambda = lambda - o.lambda;
	ans.A = A - o.A;
	ans.B = B - o.B;
	ans.omega = omega - o.omega;
	return ans;
}

double Param::h( double t) const{
	return exp(-lambda*t) *( A*cos(omega*t) + B*sin(omega*t));
}

void Param::dh(double t, double*v) const {
	v[0] = -t*exp(-lambda*t) *( A*cos(omega*t) + B*sin(omega*t));
	v[1] = exp(-lambda*t)*cos(omega*t);
	v[2] = exp(-lambda*t)*sin(omega*t);
	v[3] = exp(-lambda*t)*( -A*t*sin(omega*t) + B*t*cos(omega*t) ) ;
}

void Param::ddh(double t, double**v) const{
	double yhat = this->h(t);
	v[0][0] = t*t*exp(-lambda*t)*( A*cos(omega*t) + B*sin(omega*t) );
	v[0][1] = v[1][0] = -t*exp(-lambda*t)*cos(omega*t);
	v[0][2] = v[2][0] = -t*exp(-lambda*t)*sin(omega*t);
	v[0][3] = v[3][0] = -t*t*exp(-lambda*t)*( -A*sin(omega*t) + B*cos(omega*t) );
	v[1][1] = 0; 
	v[1][2] = v[2][1] = 0; 
	v[1][3] = v[3][1] = -t*exp(-lambda*t)*sin(omega*t);
	v[2][2] = 0;
	v[2][3] = v[3][2] = t*exp(-lambda*t)*cos(omega*t);
	v[3][3] = -t*t*exp(-lambda*t)*(A*cos(omega*t) + B*sin(omega*t));
}

void Param::disp() {
	double xi = lambda/sqrt( lambda*lambda + omega*omega);
	double w = omega/sqrt(1-xi*xi);
	cout<<"freq = "<< w/2/pi << endl
		<<"xi = "<< xi <<endl
		<<"A = "<<A<<endl
		<<"B = "<<B<<endl;
}
