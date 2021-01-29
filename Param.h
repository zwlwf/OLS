#ifndef _PARAM_H
#define _PARAM_H
struct Param {
	double lambda;
	double A;
	double B;
	int m;
	Param() {
		lambda = A = B = 0.0;
		m = 3;
	}
	Param(double l, double a, double b) {
		lambda = l; A = a; B = b; 
		m = 3;
	}

	void operator=(const Param& o);
	void operator+=(const Param& o);
	Param operator*(const double alpha);
	Param operator+(const Param& o);
	Param operator-(const Param& o);
	// args are extra parameters
	double h(double t, void* args) const;
	void dh(double t, double*v, void* args) const;
	void ddh(double t, double**v, void* args) const;
	double * vec() { // a little trick here
		return &lambda;
	}
	const double *vec() const { // a little trick here
		return &lambda;
	}
	void disp();
};

#endif
