#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <thread>
#include <iomanip>
#include "functions.h"
//#include <armadillo>
#include "constants.h"

using namespace std;
//using namespace arma;

double envelope(double x, double y, double z, double t) {
	return exp(-pow((t + z - zpulseinit), 2.0) / (2.0*sigma*sigma))+1.0*exp(-pow((t + z - zpulseinit-2.0*sigma), 2.0) / (0.5*sigma*sigma));
	//return 1.0;
}

double Fx(double x, double y, double z, double t) {

	return Estr * pow(cos(w0*(t + z)),1.0);
	//return Estr * cos(w0*(t + z - zpulseinit));
}

double Fy(double x, double y, double z, double t) {

	return poldegree*Estr * pow(cos(w0*(t + z)+phase), 1.0);
	//return Estr * cos(w0*(t + z - zpulseinit));
}

double ExDTplane(double x, double y, double z, double* rpass) {
	if (fabs(x) < *rpass) *rpass = fabs(x);

	double out = 0.0;
	for (int i = 0; i < 4; i++) {
		out += 2.0*pow(pi, 0.5)*e*a0*numberdensity*dp*a[i] / pow(B[i] + pow(u, 2.0), 0.5)*exp(-pow(x, 2.0) / (B[i] + pow(u, 2.0)))*2.0*x / (B[i] + pow(u, 2.0));
	}
	return out;
};

double dExdxDTplane(double x, double y, double z) {

	double out = 0.0;
	for (int i = 0; i < 4; i++) {
		out += 4.0*pow(pi, 0.5)*e*a0*numberdensity*dp*a[i] / pow(B[i] + pow(u, 2.0), 1.5)*exp(-pow(x, 2.0) / (B[i] + pow(u, 2.0)))*(1.0 - 2.0*pow(x, 2.0) / (B[i] + pow(u, 2.0)));
	}
	return out;
};

double UDTplane(double x, double y, double z, double* rpass) {
	if (abs(x) < *rpass) *rpass = abs(x);
	double out = 0.0;
	for (int i = 0; i < 4; i++) {
		out += 2.0*pow(pi, 0.5)*pow(e, 2.0)*a0*numberdensity*dp*a[i] / pow(B[i] + pow(u, 2.0), 0.5)*exp(-pow(x, 2.0) / (B[i] + pow(u, 2.0)));
	}
	return out;
};

double UDTplanesum(double x, double y, double z, double* rpass) {
	x = mymod(x + dp / 2.0, dp) - dp / 2.0;
	double out = 0.0;
	for (int i = -2; i < 3; i++) {
		out += UDTplane(x - i * dp, 0.0, 0.0, rpass);
	}
	return out;
}

double UDT(double x) {
	double *rpass = new double;
	*rpass = 10.0;
	double out = UDTplanesum(x, 0.0, 0.0, rpass) - UDTplanesum(dp / 2.0, 0.0, 0.0, rpass);
	delete rpass;
	return out;
}

double ExDTplanesum(double x, double y, double z, double* rpass) {
	x = mymod(x + dp / 2.0, dp) - dp / 2.0;
	double out = 0.0;
	for (int i = -2; i < 3; i++) {
		out += ExDTplane(x - i * dp, 0.0, 0.0, rpass);
	}
	return out;
}

double dExdxDTplanesum(double x, double y, double z) {
	x = mymod(x + dp / 2.0, dp) - dp / 2.0;
	double out = 0.0;
	for (int i = -2; i < 3; i++) {
		out += dExdxDTplane(x - i * dp, 0.0, 0.0);
	}
	return out;
}


double Ex(double x, double y, double z, double t) {

	// Monochromatic
	//return Estr * cos(w0*(t + z));	
	
	
	
	//Pulse
	//double phi = w0 * (t + z);	
	double phi = w0 * (t + sin(laserangle)*y + cos(laserangle)*z);
	double f = pow(sin(phi / (2.0*N)), 4.0);
	double fp = 2.0*pow(sin(phi / (2.0*N)), 3.0)*cos(phi / (2.0*N)) / N;
	/*if (phi>=0.0 && phi<=2.0*pi*N) return Estr * (f*sin(phi) - fp * cos(phi));
	else return 0.0;*/	
	return 0.0;
	
}

double Ey(double x, double y, double z, double t) {
	
	// Monochromatic
	//double phi = w0 * (t + z) - 2.0*2.0*pi+phase;
	//return poldegree*Estr * sinh(phi) / (cosh(phi)*cosh(phi));
	//return 0.0*Estr * sin(w0*(t + z));

	//Pulse
	double phi = w0 * (t + sin(laserangle)*y + cos(laserangle)*z);
	double f = pow(sin(phi / (2.0*N)), 4.0);
	double fp = 2.0*pow(sin(phi / (2.0*N)), 3.0)*cos(phi / (2.0*N)) / N;
	if (phi>=0.0 && phi<=2.0*pi*N) return -cos(laserangle)*Estr * (f*sin(phi) - fp * cos(phi));
	else return 0.0;
	//return 0.0;
}

double Ez(double x, double y, double z, double t) {
	double phi = w0 * (t + sin(laserangle)*y + cos(laserangle)*z);
	double f = pow(sin(phi / (2.0*N)), 4.0);
	double fp = 2.0*pow(sin(phi / (2.0*N)), 3.0)*cos(phi / (2.0*N)) / N;
	if (phi >= 0.0 && phi <= 2.0*pi*N) return sin(laserangle)*Estr * (f*sin(phi) - fp * cos(phi));
	else return 0.0;
}

double dExdx(double x, double y, double z) {
	//return -(phifinal(x+potfine/2.0,y,z)-phifinal(x-potfine/2.0,y,z))/potfine;
	//return dExdxDTstringsum(x,y,z);
	return dExdxDTplanesum(x - z * z / (2.0*R), y, z);;
}

double dExdy(double x, double y, double z) {
	//return -(phifinal(x+potfine/2.0,y,z)-phifinal(x-potfine/2.0,y,z))/potfine;
	//return dExdyDTstringsum(x,y,z);
	return 0.0;
}

double dEydy(double x, double y, double z) {
	//return -(phifinal(x+potfine/2.0,y,z)-phifinal(x-potfine/2.0,y,z))/potfine;
	//return dEydyDTstringsum(x,y,z);
	return 0.0;
}

double Bx(double x, double y, double z, double t) {
	return -sin(laserangle)*Ez(x, y, z, t)+cos(laserangle)*Ey(x, y, z, t);
}

double By(double x, double y, double z, double t) {
	return -cos(laserangle)*Ex(x, y, z, t);
}

double Bz(double x, double y, double z, double t) {
	return sin(laserangle)*Ex(x, y, z, t);
}
double Ex2(double x, double y, double z, double t, double Eparam) {
	return Eparam *( cos(w0*(t + z)));
}

double Ey2(double x, double y, double z, double t, double Eparam) {
	//return -(phifinal(x,y+potfine/2.0,z)-phifinal(x,y-potfine/2.0,z))/potfine;
	//return EyDTstringsum(x,y,z);
	return Eparam * pow(cos(w0*(t + z)+phase), 1.0);
}

double Ez2(double x, double y, double z, double t, double Eparam) {
	//return -(phifinal(x,y,z+potfine/2.0)-phifinal(x,y,z-potfine/2.0))/potfine;
	return 0.0;
}


double Bx2(double x, double y, double z, double t, double Eparam) {
	return Ey2(x, y, z, t, Eparam);
}

double By2(double x, double y, double z, double t, double Eparam) {
	return -Ex2(x, y, z, t, Eparam);
}

double Bz2(double x, double y, double z, double t, double Eparam) {
	return 0.0;
}

int func(double t, const double y[], double f[],
	void *params)
{
	double ptemp;
	double gammatemp;
	double *par = (double *)params;
	double gam = par[0];
	double *rpass = par + 1;

	//ptemp=pow(pow(y[0],2)+pow(y[1],2)+pow(y[2]+E*v0,2),0.5);
	gammatemp = y[6] + gam;
	double v0 = 1 - 0.5*pow(gam, -2.0);

	double crystalstat = 1.0;

	//double Fx=q*(Ex(t,y[3],y[4],y[5]+v0*t,y[7])-(y[2]+v0)*By(t,y[3],y[4],y[5]+v0*t,y[7]));
	//double Fy=q*(Ey(t,y[3],y[4],y[5]+v0*t,y[7])+(y[2]+v0)*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));
	//double Fz=q*(y[0]*By(t,y[3],y[4],y[5]+v0*t,y[7])-y[1]*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));

	double Fx = q * (Ex(y[3], y[4], y[5] + v0 * t, t) + y[1]*Bz(y[3], y[4], y[5] + v0 * t, t) - (y[2] + v0)*By(y[3], y[4], y[5] + v0 * t,t));
	double Fy = q * (Ey(y[3], y[4], y[5] + v0 * t,t) + (y[2] + v0)*Bx(y[3], y[4], y[5] + v0 * t,t)- y[0] * Bz(y[3], y[4], y[5] + v0 * t, t));
	double Fz = q * Ez(y[3], y[4], y[5] + v0 * t,t) + q * (y[0] * By(y[3], y[4], y[5] + v0 * t,t) - y[1] * Bx(y[3], y[4], y[5] + v0 * t,t));

	
	//ADDING RADIATION REACTION
	/*if (rrstat == 1) {

		double chicur = pow(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0), 0.5)*gammatemp*e / (m*m);
		double baierfac = 1.0 / pow(1.0 + 4.8*(1.0 + chicur)*log(1.0 + 1.7*chicur) + 2.44*pow(chicur, 2.0), 2.0 / 3.0);
		baierfac = 1.0;
		Fx += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[0] * dExdx(y[3], y[4], y[5] + v0 * t) + y[1] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ex(y[3], y[4], y[5] + v0 * t, rpass)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[0] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));

		Fy += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[1] * dEydy(y[3], y[4], y[5] + v0 * t) + y[0] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ey(y[3], y[4], y[5] + v0 * t)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[1] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));

		Fz += baierfac * (-2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*(y[2] + v0)*(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
	}*/
	// DONE RR

	double gamprime = (y[0] * Fx + y[1] * Fy + (y[2] + v0)*Fz) / m;
	//gamprime=0.0;

	f[0] = Fx / (gammatemp*m) - gamprime / gammatemp * y[0];
	f[1] = Fy / (gammatemp*m) - gamprime / gammatemp * y[1];	
	f[2] = 1.0 / (gammatemp*m)*(Fz*(y[0] * y[0] + y[1] * y[1] + 1.0 / (gammatemp*gammatemp)) - Fx * y[0] - Fy * y[1]);

	f[3] = y[0];
	f[4] = y[1];	
	f[5] = 0.5*(pow(gam, -2.0) - pow(gammatemp, -2.0) - pow(y[0], 2.0) - pow(y[1], 2.0));

	


	f[6] = gamprime;


	return GSL_SUCCESS;
}

//int func2(double t, const double y[], double f[],
//	void *params)
//{
//	double ptemp;
//	double gammatemp;
//	double *par = (double *)params;
//	double gam = par[0];
//	double *rpass = par + 1;
//	double Eparamx = par[1];
//	double Eparamy = par[2];
//
//	//ptemp=pow(pow(y[0],2)+pow(y[1],2)+pow(y[2]+E*v0,2),0.5);
//	gammatemp = y[6] + gam;
//	double v0 = 1 - 0.5*pow(gam, -2.0);
//
//
//
//	//double Fx=q*(Ex(t,y[3],y[4],y[5]+v0*t,y[7])-(y[2]+v0)*By(t,y[3],y[4],y[5]+v0*t,y[7]));
//	//double Fy=q*(Ey(t,y[3],y[4],y[5]+v0*t,y[7])+(y[2]+v0)*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));
//	//double Fz=q*(y[0]*By(t,y[3],y[4],y[5]+v0*t,y[7])-y[1]*Bx(t,y[3],y[4],y[5]+v0*t,y[7]));
//
//	double Fx = q * (Ex2(y[3], y[4], y[5] + v0 * t, t, Eparamx) - (y[2] + v0)*By2(y[3], y[4], y[5] + v0 * t, t, Eparamx));
//	double Fy = q * (Ey2(y[3], y[4], y[5] + v0 * t, t, Eparamy) + (y[2] + v0)*Bx2(y[3], y[4], y[5] + v0 * t, t, Eparamy));
//	double Fz = q * Ez2(y[3], y[4], y[5] + v0 * t, t, Eparamx) + q * (y[0] * By2(y[3], y[4], y[5] + v0 * t, t, Eparamx) - y[1] * Bx2(y[3], y[4], y[5] + v0 * t, t, Eparamy));
//
//
//	//ADDING RADIATION REACTION
//	//if (rrstat == 1) {
//	//	double chicur = pow(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0), 0.5)*gammatemp*e / (m*m);
//	//	double baierfac = 1.0 / pow(1.0 + 4.8*(1.0 + chicur)*log(1.0 + 1.7*chicur) + 2.44*pow(chicur, 2.0), 2.0 / 3.0);
//	//	//baierfac=1.0;
//	//	Fx += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[0] * dExdx(y[3], y[4], y[5] + v0 * t) + y[1] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ex(y[3], y[4], y[5] + v0 * t, rpass)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[0] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
//
//	//	Fy += baierfac * (1.0*2.0*pow(q, 3.0) / (3.0*m)*gammatemp*(y[1] * dEydy(y[3], y[4], y[5] + v0 * t) + y[0] * dExdy(y[3], y[4], y[5] + v0 * t)) + 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*Ey(y[3], y[4], y[5] + v0 * t)*(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t)) - 2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*y[1] * (pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
//
//	//	Fz += baierfac * (-2.0 / 3.0*pow(q, 4.0) / pow(m, 2.0)*pow(gammatemp, 2.0)*(y[2] + v0)*(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0) - pow(y[0] * Ex(y[3], y[4], y[5] + v0 * t, rpass) + y[1] * Ey(y[3], y[4], y[5] + v0 * t), 2.0)));
//	//}
//	// DONE RR
//
//	double gamprime = (y[0] * Fx + y[1] * Fy + (y[2] + v0)*Fz) / m;
//	//gamprime=0.0;
//
//	f[0] = Fx / (gammatemp*m) - gamprime / gammatemp * y[0];
//	f[1] = Fy / (gammatemp*m) - gamprime / gammatemp * y[1];
//	f[2] = Fz / (gammatemp*m) - gamprime / gammatemp * (y[2] + v0);
//
//	f[3] = y[0];
//	f[4] = y[1];
//	f[5] = y[2];
//	f[6] = gamprime;
//
//
//	return GSL_SUCCESS;
//}

void trajecsolver(double* upup1, double* downdown1, double* updown1, double* downup1, double* upup2, double* downdown2, double* updown2, double* downup2, double Em,  int curN, int Eindex, double thetax_length, double thetay_length) {

	
	
	double gam0 = Em / m;
	//introduce the photon to be converted
	double *k = new double[4];
	k[0] = photE;
	k[1] = thetax_phot * k[0];
	k[2] = thetay_phot * k[0];
	k[3] = (1.0 - 0.5*(thetax_phot*thetax_phot + thetay_phot * thetay_phot))*k[0];

	double *hvec = new double[4];
	double *polphot1 = new double[4]();
	double *polphot2 = new double[4]();
	double thetapol = 0.0;
	hvec[0] = 0.0;
	hvec[1] = sin(thetapol);
	hvec[2] = -cos(thetapol);
	hvec[3] = 0.0;
	mycross(k, hvec, polphot1);
	normalize(polphot1);

	thetapol = pi/2.0;
	hvec[0] = 0.0;
	hvec[1] = sin(thetapol);
	hvec[2] = -cos(thetapol);
	hvec[3] = 0.0;
	mycross(k, hvec, polphot2);
	normalize(polphot2);


	double nx, ny, nz;
	nx = thetax_phot;
	ny = thetay_phot;
	nz = (1.0 - 0.5*(thetax_phot*thetax_phot + thetay_phot * thetay_phot));
	double *nvec = new double[4];
	nvec[0] = 0.0;
	nvec[1] = nx;
	nvec[2] = ny;
	nvec[3] = nz;

	double *ecrossn1 = new double[4]();
	double *ecrossn2 = new double[4]();
	ecrossn1[0] = 0.0;
	ecrossn2[0] = 0.0;
	mycross(polphot1, nvec, ecrossn1);
	mycross(polphot2, nvec, ecrossn2);


	complex<double> sm_u;
	complex<double> sm_d;
	complex<double> sp_u;
	complex<double> sp_d;

	for (int j = 0; j < loopsize; j++) {

		double *trajec = new double[8 * pointsprim]();

		int index = j + curN * loopsize;
		int index_x = index % thetafine;
		int index_y = (index - index_x) / thetafine;



		ofstream orbit;
		char a[150];
		char b[150];
		strcpy_s(a, "orbit");
		sprintf_s(b, "%i", Eindex);	
		strcat_s(a, b);
		sprintf_s(b, "-%i.txt", index);
		strcat_s(a, b);

		if (writestat == 1) orbit.open(a);
		if (writestat == 1) orbit << scientific << setprecision(16);

		const gsl_rng_type * Trand;
		gsl_rng * r;
		gsl_rng_env_setup();

		Trand = gsl_rng_default;
		r = gsl_rng_alloc(Trand);
		gsl_rng_set(r, curN);

		double t0 = 0.0;
		double *params = new double[2];
		params[0] = Em / m;
		double *rpass = params + 1;
		double t = t0;
		int points = pointsprim;
		
		//double *out = new double[freqiters]();


		double vx = thetax_phot - 0.5*thetax_length + (double)index_x / ((double)thetafine - 1.0)*thetax_length;
		double vy = thetay_phot - 0.5*thetay_length + (double)index_y / ((double)thetafine - 1.0)*thetay_length;
		double xpos = 0.0;
		double ypos = 0.0;

		double y[7];
		y[0] = vx; //vx
		y[1] = vy; //vy
		y[2] = (-pow(vx, 2.0) - pow(vy, 2.0)) / 2.0;
		y[3] = xpos;
		y[4] = ypos;
		y[5] = 0.0;
		y[6] = 0.0;
		double gam = Em / m;

		double gamtemp = y[6] + gam;
		trajec[0 * points] = t;
		trajec[1 * points] = y[0];
		trajec[2 * points] = y[1];
		trajec[3 * points] = y[2];
		trajec[4 * points] = y[3];
		trajec[5 * points] = y[4];
		trajec[6 * points] = y[5];
		trajec[7 * points] = y[6];

		double thetax_max = y[0];
		double thetax_min = y[0];
		double thetay_max = y[1];
		double thetay_min = y[1];

		gsl_odeiv2_system sys = { func, NULL, 7, params };
		gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-8, 1e-8, 1e-8);

		complex<double> Ix(0.0, 0.0);
		complex<double> Iy(0.0, 0.0);
		complex<double> J(0.0, 0.0);
		complex<double> K1_x(0.0, 0.0);
		complex<double> K1_y(0.0, 0.0);
		complex<double> K1_z(0.0, 0.0);
		complex<double> K2_x(0.0, 0.0);
		complex<double> K2_y(0.0, 0.0);
		complex<double> K2_z(0.0, 0.0);

		int i, status;
		double ti;
		double dt = tprim / ((double)points - 1);
		complex<double> dtc(dt, 0.0);
		for (i = 1; i <= points - 1; i++)
		{
			*rpass = aL;
			ti = t0 + (double)i * dt;
			while (t < ti) {
				//status = gsl_odeiv_evolve_apply (e, c, s,&sys,&t,ti,&h, y);
				status = gsl_odeiv2_driver_apply(d, &t, ti, y);
				if (status != GSL_SUCCESS)
				{
					printf("error, return value=%d\n", status);
					break;
				}
			}
			//   printf ("%.5e %.5e %.5e\n", t, y[3], y[4], y[5]);


			//double scatvar = fmod(abs(y[3]), dp);
			//double Gauss = dp / (u*pow(2.0*pi, 0.5))*(exp(-pow(scatvar, 2.0) / (2.0*u*u)) + exp(-pow(dp - scatvar, 2.0) / (2.0*u*u)));
			//double nel = -1.0*dExdx(scatvar, 0.0, 0.0) / (4.0*pi*pow(alfa, 0.5)) + Z * numberdensity*Gauss;
			//double sigsqncl = Gauss * pow(cscat / (m*(gam + y[6])), 2.0)*l_f / (double)pointssec;
			//double sigsqel = 1.0*pi*pow(alfa, 2.0) / pow(m*(gam + y[6]), 2.0)*(log(2.0*m*pow(gam + y[6], 2.0) / ISi) - 1.0) * nel *l_f / (double)pointssec;
			////sigsqncl = pow(cscat / (m*(gam + y[6])), 2.0)*tprim / (double)pointssec;
			////sigsqel = 1.0*pi*pow(alfa, 2.0) / pow(m*(gam + y[6]), 2.0)*(log(2.0*m*pow(gam + y[6], 2.0) / ISi) - 1.0) * Z*numberdensity*tprim / (double)pointssec;
			////cout << sigsqncl / sigsqel << endl;
			//double scatx = pow(1.0*sigsqel + 1.0*sigsqncl, 0.5);
			//double scaty = scatx;

			//scatx = 0.0*scatx*gsl_ran_gaussian(r, 1.0);
			//scaty = 0.0*scaty*gsl_ran_gaussian(r, 1.0);
			////y[2] = y[2] - 0.5*(pow(scatx, 2.0) + pow(scaty, 2.0)) - y[0] * scatx - y[1] * scaty;
			//y[0] = y[0] + scatx;
			//y[1] = y[1] + scaty;
			//y[2] = 0.5*(pow(gam, -2.0) - pow(gamtemp, -2.0) - pow(y[0], 2.0) - pow(y[1], 2.0));

			//if (0) {

			//	double gamold = y[6] + gam;
			//	//cout << scatx * gamold << endl;
			//	double deltagam = -2.0*alfa / (3.0*pi)*pow(gamold, 3.0)*(pow(scatx, 2.0) + pow(scaty, 2.0))*3.0 / 4.0;
			//	double vrat = 1.0 + deltagam / pow(gamold, 3.0);
			//	y[0] = y[0] * vrat;
			//	y[1] = y[1] * vrat;
			//	y[2] += deltagam / pow(gamold, 3.0);
			//	y[6] += deltagam;
			//}

			//double field = pow(pow(Ex(y[3], y[4], y[5] + v0 * t, rpass), 2.0) + pow(Ey(y[3], y[4], y[5] + v0 * t), 2.0), 0.5);
			//double chi = gamtemp * field*e / pow(m, 2.0);

			//double Emek = UDT(y[3]) + 0.5*gamtemp*m*y[0] * y[0];

			if (writestat == 1 && i % 1 == 0) orbit << y[3] << " " << y[4] << " " << y[5] << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[6] + gam << " " << t << endl;

			trajec[i + 0 * points] = t;
			trajec[i + 1 * points] = y[0];
			trajec[i + 2 * points] = y[1];
			trajec[i + 3 * points] = y[2];
			trajec[i + 4 * points] = y[3];
			trajec[i + 5 * points] = y[4];
			trajec[i + 6 * points] = y[5];
			trajec[i + 7 * points] = y[6];

			double x = trajec[4 * points + i];
			double y = trajec[5 * points + i];
			double z = trajec[6 * points + i];

			vx = trajec[1 * points + i];
			vy = trajec[2 * points + i];
			double vz = trajec[3 * points + i];
			
			double betaxdot, betaydot, betazdot;

			
			betaxdot = (trajec[1 * points + i] - trajec[1 * points + i-1]) / dt;
			betaydot = (trajec[2 * points + i] - trajec[2 * points + i-1]) / dt;
			betazdot = (trajec[3 * points + i] - trajec[3 * points + i-1]) / dt;
			

			complex<double> exponent(0.0, photE*Em/(photE-Em)*(-nx * x - ny * y - z + 0.5*(nx*nx + ny * ny + 1.0 / (gam0*gam0))*t));
			complex<double> faktor(1.0 / ((0.5*(1.0 / (gam0*gam0) + nx * nx + ny * ny) - nx * vx - ny * vy - vz)*(0.5*(1.0 / (gam0*gam0) + nx * nx + ny * ny) - nx * vx - ny * vy - vz)), 0.0);
			
			complex<double> faktor2((0.5*(1.0 / (gam0*gam0) + nx * nx + ny * ny) - nx * vx - ny * vy - vz), 0.0);

			complex<double> v1(ny*(betaydot*(nx - vx) - betaxdot * (ny - vy)) - betaxdot * 0.5*(1.0 / (gam0*gam0) - nx * nx - ny * ny - 2.0*vz) + (nx - vx)*betazdot, 0.0);
			complex<double> v2((ny - vy)*betazdot - betaydot * 0.5*(1.0 / (gam0*gam0) - nx * nx - ny * ny - 2.0*vz) - nx * ((nx - vx)*betaydot - (ny - vy)*betaxdot), 0.0);
			v1 = v1 * faktor*exp(exponent);
			v2 = v2 * faktor*exp(exponent);
			complex<double> vJ(nx*(betaxdot)+(ny)*(betaydot)+betazdot, 0.0);
			vJ = vJ * faktor*exp(exponent);
			
			double edotv1 = polphot1[1] * vx + polphot1[2] * vy + polphot1[3];
			double edota1 = polphot1[1] * betaxdot + polphot1[2] * betaydot + polphot1[3]*betazdot;

			double edotv2 = polphot2[1] * vx + polphot2[2] * vy + polphot2[3];
			double edota2 = polphot2[1] * betaxdot + polphot2[2] * betaydot + polphot2[3] * betazdot;

			double ndota = nx * betaxdot + ny * betaydot + betazdot;

			complex<double> k1_x = (betaxdot*edotv1+vx*edota1)*faktor2+vx*edotv1*ndota;
			complex<double> k1_y = (betaydot*edotv1 + vy * edota1)*faktor2 + vy * edotv1*ndota;
			complex<double> k1_z = (betazdot*edotv1 + edota1)*faktor2 + edotv1*ndota;

			complex<double> k2_x = (betaxdot*edotv2 + vx * edota2)*faktor2 + vx * edotv2*ndota;
			complex<double> k2_y = (betaydot*edotv2 + vy * edota2)*faktor2 + vy * edotv2*ndota;
			complex<double> k2_z = (betazdot*edotv2 + edota2)*faktor2 + edotv2 * ndota;

			k1_x = k1_x * faktor*exp(exponent);
			k1_y = k1_y * faktor*exp(exponent);
			k1_z = k1_z * faktor*exp(exponent);
			
			k2_x = k2_x * faktor*exp(exponent);
			k2_y = k2_y * faktor*exp(exponent);
			k2_z = k2_z * faktor*exp(exponent);

			Ix += v1 * dt;
			Iy += v2 * dt;
			J += vJ * dt;

			K1_x += k1_x * dt;
			K1_y += k1_y * dt;
			K1_z += k1_z * dt;

			K2_x += k2_x * dt;
			K2_y += k2_y * dt;
			K2_z += k2_z * dt;
			


		}//end time step
		complex<double> imagu(0.0, 1.0);
		complex<double> nullc(0.0, 0.0);

		complex<double> A = -photE*Em * (Ix*ecrossn1[1]+ Iy * ecrossn1[2]);
		//complex<double> Idotn = Ix * nx + Iy * ny;
		complex<double> Idotn(0.0, 0.0);
		complex<double> poldotI = polphot1[1] * Ix + polphot1[2] * Iy;

		complex<double> Bx = polphot1[1] * photE*(m*J + Em * Idotn) - 1.0*2.0*Em*Em*K1_x - Em * photE*nx*poldotI;
		complex<double> By = polphot1[2] * photE*(m*J + Em * Idotn) - 1.0*2.0*Em*Em*K1_y - Em * photE*ny*poldotI;
		complex<double> Bz = polphot1[3] * photE*(m*J + Em * Idotn) - 1.0*2.0*Em*Em*K1_z - Em * photE*poldotI;


		complex<double> M11 = A - imagu * Bz;
		complex<double> M12 = -imagu * Bx - By;
		complex<double> M21 = -imagu * Bx + By;
		complex<double> M22 = A + imagu * Bz;

		//Spin along z up-up
		sm_u = complex<double>(1.0, 0.0);
		sm_d = complex<double>(0.0, 0.0);
		sp_u = complex<double>(0.0, 0.0);
		sp_d = complex<double>(1.0, 0.0);
		upup1[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));
		
		//Spin along z down-down
		sm_u = complex<double>(0.0, 0.0);
		sm_d = complex<double>(1.0, 0.0);
		sp_u = complex<double>(1.0, 0.0);
		sp_d = complex<double>(0.0, 0.0);

		downdown1[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));

		//Spin along z up-down
		sm_u = complex<double>(1.0, 0.0);
		sm_d = complex<double>(0.0, 0.0);
		sp_u = complex<double>(1.0, 0.0);
		sp_d = complex<double>(0.0, 0.0);
		updown1[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));

		//Spin along z down-up
		sm_u = complex<double>(0.0, 0.0);
		sm_d = complex<double>(1.0, 0.0);
		sp_u = complex<double>(0.0, 0.0);
		sp_d = complex<double>(1.0, 0.0);
		downup1[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));

		A = -photE * Em * (Ix*ecrossn2[1] + Iy * ecrossn2[2]);
		//complex<double> Idotn = Ix * nx + Iy * ny;
		
		poldotI = polphot2[1] * Ix + polphot2[2] * Iy;

		Bx = polphot2[1] * photE*(m*J + Em * Idotn) - 1.0*2.0*Em*Em*K2_x - Em * photE*nx*poldotI;
		By = polphot2[2] * photE*(m*J + Em * Idotn) - 1.0*2.0*Em*Em*K2_y - Em * photE*ny*poldotI;
		Bz = polphot2[3] * photE*(m*J + Em * Idotn) - 1.0*2.0*Em*Em*K2_z - Em * photE*poldotI;

		M11 = A - imagu * Bz;
		M12 = -imagu * Bx - By;
		M21 = -imagu * Bx + By;
		M22 = A + imagu * Bz;

		//Spin along z up-up
		sm_u = complex<double>(1.0, 0.0);
		sm_d = complex<double>(0.0, 0.0);
		sp_u = complex<double>(0.0, 0.0);
		sp_d = complex<double>(1.0, 0.0);
		upup2[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));

		//Spin along z down-down
		sm_u = complex<double>(0.0, 0.0);
		sm_d = complex<double>(1.0, 0.0);
		sp_u = complex<double>(1.0, 0.0);
		sp_d = complex<double>(0.0, 0.0);

		downdown2[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));

		//Spin along z up-down
		sm_u = complex<double>(1.0, 0.0);
		sm_d = complex<double>(0.0, 0.0);
		sp_u = complex<double>(1.0, 0.0);
		sp_d = complex<double>(0.0, 0.0);
		updown2[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));

		//Spin along z down-up
		sm_u = complex<double>(0.0, 0.0);
		sm_d = complex<double>(1.0, 0.0);
		sp_u = complex<double>(0.0, 0.0);
		sp_d = complex<double>(1.0, 0.0);
		downup2[index] = norm(conj(sm_u)*(M11*sp_u + M12 * sp_d) + conj(sm_d)*(M21*sp_u + M22 * sp_d));


		if (writestat == 1) orbit.close();



		gsl_odeiv2_driver_free(d);
		gsl_rng_free(r);

		delete[] params;
		delete[] trajec;

	}

}



int myintmod(int a, int b) {
	int x;
	if (a > 0) x = a % b;
	else if (a < 0) x = a % b + b;
	else x = 0;
	return x;
};



double lnfac(int n) {
double out;
double x=(double)n;

if (n!=0) {
out=x*log(x)-x+0.5*log(2.0*pi*x)+1.0/(12.0*x)-1.0/( 360.0*pow(x,3.0) );
}
else out=0;

return out;
};

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
};

double mylag(int n, int alfa, double x) {
double out=0.0;
for (int i=0; i<=n; i++) {
out+=pow(-1.0,i)*factorial(n+alfa)/factorial(alfa+i)/factorial(n-i)*pow(x,(double)i)/factorial(i);
};

return out;
}

double trapz (double *t, double *f, int points) {
double dummy=0;
for (int i=0; i<points-1; i++) {
if ( (t[i+1]-t[i])>0.0 ) dummy+=(t[i+1]-t[i])*(f[i+1]+f[i])/2.0;
}
return dummy;
};

void mycross (double *a, double *b, double *c) {
c[1]=a[2]*b[3]-a[3]*b[2];
c[2]=a[3]*b[1]-a[1]*b[3];
c[3]=a[1]*b[2]-a[2]*b[1];
};

void normalize (double *a) {
double x=pow(a[1]*a[1]+a[2]*a[2]+a[3]*a[3],0.5);
a[1]=a[1]/x;
a[2]=a[2]/x;
a[3]=a[3]/x;
};

double mymod(double a, double b) {
	double x;
	if (a >= 0) x = fmod(a, b);
	else x = fmod(a, b) + b;
	return x;
};



int sgn(double a) {
	if (a >= 0) return 1;
	else return -1;
}

double fdot(double *x1, double *x2) {
double out=x1[0]*x2[0];
for (int i=1; i<4; i++) {
out-=x1[i]*x2[i];
}
return out;
}






double myinterp(double *a, double *b, double x, int points) {


	bool success = 0;

	if (x > a[points - 1]) {
		return b[points - 1];
	}

	if (x < a[0]) {
		return b[0];
	}

	else {
		for (int i = 0; i < points - 1; i++) {
			if (x >= a[i] && x < a[i + 1]) {
				success = 1;
				return b[i] + (b[i + 1] - b[i]) / (a[i + 1] - a[i])*(x - a[i]);
			}
		}
	}

	if (success == 0) {
		return 0;
	}
}