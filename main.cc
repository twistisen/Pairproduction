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
#include "constants.h"
#include <random>

using namespace std;


int main() {
	
	char a[150];
	char b[150];
	ofstream enh1,enh2,enh3,file;
	strcpy_s(a, "spectrum.txt");
	enh1.open(a);
	
	int div = 500;
	file.open("potential.txt");
	for (int i = 0; i < div; i++) {
		double x = (double)i / (double)div*dp;
		file << x << " " << UDT(x) << endl;
	}
	file.close();

	
	double *resultupup1 = new double[Efine]();
	double *resultdowndown1 = new double[Efine]();
	double *resultupdown1 = new double[Efine]();
	double *resultdownup1 = new double[Efine]();

	double *resultupup2 = new double[Efine]();
	double *resultdowndown2 = new double[Efine]();
	double *resultupdown2 = new double[Efine]();
	double *resultdownup2 = new double[Efine]();
	

	for (int l = 0; l < Efine; l++) {

		

		cout << l << endl;
		double Em = Emin + (double)l / (double)Efine*(Emax - Emin);

		const double thetax_length = m / Em * (1.5 + 2.0*eta);
		const double thetay_length = m / Em * (2.5 + 0.0*eta);
		
		
		
		//electron first, positron second
		double *upup1 = new double[thetafine*thetafine]();
		double *downdown1 = new double[thetafine*thetafine]();
		double *updown1 = new double[thetafine*thetafine]();
		double *downup1 = new double[thetafine*thetafine]();

		double *upup2 = new double[thetafine*thetafine]();
		double *downdown2 = new double[thetafine*thetafine]();
		double *updown2 = new double[thetafine*thetafine]();
		double *downup2 = new double[thetafine*thetafine]();
		
		char a[150];
		char b[150];
		strcpy_s(a, "angle1-");
		sprintf_s(b, "%i.txt", l);
		strcat_s(a, b);
		if (writeangle == 1) enh2.open(a);

		strcpy_s(a, "angle2-");
		sprintf_s(b, "%i.txt", l);
		strcat_s(a, b);
		if (writeangle == 1) enh3.open(a);

		int multiplum = thetafine * thetafine / threadN;
		std::vector < std::thread > threads(threadN);
		
			for (int i = 0; i < threadN; i++) {
				threads[i] = std::thread(trajecsolver,upup1,downdown1,updown1,downup1, upup2, downdown2, updown2, downup2,Em, i, l, thetax_length,thetay_length);
			}

			for (int i = 0; i < threadN; i++) {
				threads[i].join();
			}

			for (int i = 0; i < thetafine; i++) {
				for (int j = 0; j < thetafine; j++) {
					int index = j + i * thetafine;
					int index_x = index % thetafine;
					int index_y = (index - index_x) / thetafine;
					double upupcontrib = upup1[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					double downdowncontrib = downdown1[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					double updowncontrib = updown1[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					double downupcontrib = downup1[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					resultupup1[l] += upupcontrib;
					resultdowndown1[l] += downdowncontrib;
					resultupdown1[l] += updowncontrib;
					resultdownup1[l] += downupcontrib;
					if (writeangle == 1) enh2 << upupcontrib + downdowncontrib + updowncontrib + downupcontrib << " ";

					upupcontrib = upup2[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					downdowncontrib = downdown2[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					updowncontrib = updown2[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					downupcontrib = downup2[index] * e*e / (4.0*pi*pi*photE*photE*photE*Em*Em*4.0)*thetax_length / ((double)thetafine - 1.0)*thetay_length / ((double)thetafine - 1.0);
					resultupup2[l] += upupcontrib;
					resultdowndown2[l] += downdowncontrib;
					resultupdown2[l] += updowncontrib;
					resultdownup2[l] += downupcontrib;
					if (writeangle == 1) enh3 << upupcontrib + downdowncontrib + updowncontrib + downupcontrib << " ";

				}
				if (writeangle == 1) enh2 << endl;
				if (writeangle == 1) enh3 << endl;
			}

			if (writeangle == 1) enh2.close();
			if (writeangle == 1) enh3.close();

			enh1 << Em << " " << resultupup1[l] << " " << resultdowndown1[l] << " " << resultupdown1[l] << " " << resultdownup1[l] << " " << resultupup2[l] << " " << resultdowndown2[l] << " " << resultupdown2[l] << " " << resultdownup2[l] << endl;

			
		
		delete[] upup1;
		delete[] downdown1;
		delete[] updown1;
		delete[] downup1;
		delete[] upup2;
		delete[] downdown2;
		delete[] updown2;
		delete[] downup2;
	}
		
	enh1.close();


}