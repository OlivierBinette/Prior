#ifndef SRC_RNG_GAMMAVAR_H_
#define SRC_RNG_GAMMAVAR_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "../utils.h"
#include "NormalVar.h"
#include <random>

using namespace std;

/*
struct GammaVar {
	double alpha, beta;
	gamma_distribution<double> dist;
	mt19937 gen;

	GammaVar(double alph, double bet) {
		alpha = alph; beta = bet;
		dist = gamma_distribution<double>(alpha, beta);
		gen = mt19937();
	}

	GammaVar(const GammaVar& o)  {
		alpha = o.alpha;
		beta = o.beta;
		dist = gamma_distribution<double>(alpha, beta);
		gen = mt19937();
	}

	double var() {
		return dist(gen);
	}
};*/

/**
 * Structure for gamma variates (Numerical Recipes, 3rd ed., p.370).
 *
 * USAGE
 * =====
 * >>> double alpha = 0.2, beta = 1.3;
 * >>> Gammavar gam(alpha, beta);   // Generator of Gamma(alpha, beta) variates.
 * >>> gam.var();			  		// Returns a Gamma variate.
 *
 */
struct GammaVar : NormalVar {
	double alpha, oalpha, beta;
	double a1, a2;

	/* Structure for gamma deviates. */
	GammaVar(double alph, double bet)
	: NormalVar(0., 1.), alpha(alph), oalpha(alph), beta(bet) {
		if (alpha <= 0.) throw invalid_argument("Non-positive alpha in GammaVar.");
		if (alpha < 1.) alpha += 1.;
		a1 = alpha - 1./3.;
		a2 = 1./sqrt(9.*a1);
	}

	GammaVar(const GammaVar& o) : NormalVar(o) {
		alpha = o.alpha; oalpha = o.oalpha; beta = o.beta;
		a1 = o.a1; a2 = o.a2;
	}

	/** Returns a gamma deviate. */
	double var() {
		double u,v,x;
		do {
			do {
				x = NormalVar::var();
				v = 1. + a2*x;
			} while (v <= 0.);
			v = v*v*v;
			u = doub();
		} while(u > 1. - 0.331*x*x*x*x && log(u) > 0.5*x*x + a1*(1. - v + log(v)));
		if (alpha == oalpha) return a1*v/beta;
		else {
			do u = doub(); while (u == 0.);
			return pow(u, 1./oalpha)*a1*v/beta;
		}
	}
};

/**
 * Informal test of <GammaVar>. Prints readable test cases and results in the console.
 */
void test_GammaVar() {
	cout << "Testing Gammadev." << endl << "-----------------" << endl;

	double alpha = 1.2, beta = 0.7;
	GammaVar gamma(alpha, beta);

	// Filling a vector <a> with <N> random entries.
	int N = 10000;
	double* a = new double[N];
	for(int i = 0; i < N; i++) a[i] = gamma.var();

	// Calculating the empirical mean.
	double s = sum(a, N);
	s /= N;
	// Printing the result.
	cout << "Using the first " << N << " elements of the Gammadev("<<alpha<<", "<<beta<<") sequence, we have the"<<endl
		<< " empirical mean " << s << " versus the expected mean " << alpha/beta <<" . The seed is " << Ran::seed << " ." <<endl<<endl;

	// Testing other parameters.
	alpha = 0.7, beta = 1.3;
	gamma = GammaVar(alpha, beta);
	for(int i = 0; i < N; i++) a[i] = gamma.var();
	s = sum(a, N);
	s /= N;
	// Printing the result.
	cout << "Using the first " << N << " elements of the Gammadev("<<alpha<<", "<<beta<<") sequence, we have the"<<endl
		<< " empirical mean " << s << " versus the expected mean " << alpha/beta <<" . The seed is " << Ran::seed << " ."<<endl<<endl;


	// Calculating the (uncorrected) empirical variance.
	double v = 0.0;
	for (int i = 0; i < N; i++) v += (a[i] - s)*(a[i] - s);
	v /= N;
	//Printing the results.
	cout << "Using the first " << N << " elements of the Gammadev("<<alpha<<", "<<beta<<") sequence, we have the"<<endl
		<< " empirical variance " << v << " versus the expected variance " << alpha/(beta*beta) <<" . The seed is " << Ran::seed << " ."<<endl<<endl;
}



#endif /* SRC_RNG_GAMMAVAR_H_ */
