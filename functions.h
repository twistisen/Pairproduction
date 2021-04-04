#ifndef FUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FUNCTIONS_H
#include <complex>


using namespace std;

double envelope(double, double, double, double);

double lnfac(double);
int factorial(int);
double trapz (double*, double*, int);
double mymod(double, double);
int myintmod(int, int);
int sgn(double);
void mycross(double*,double*,double*);
void normalize(double*);
double myinterp(double*, double*, double, int);

double fdot(double*, double*);
void trajecsolver(double*, double*, double*, double*, double*, double*, double*, double*, double, int, int, double, double);
int func(double, const double y[], double f[],void*);
int func2(double, const double y[], double f[], void*);

double Ex(double, double, double, double*);
double Ey(double, double, double);
double Ez(double, double, double);
double UDT(double);
double dExdx(double, double, double);


//thrust::complex<double> M_dev(double*, double*, int, int, double, double, int, double*, double*, double*, double*, int, int, int);


#endif
