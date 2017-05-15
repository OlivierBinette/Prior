#ifndef SRC_RNG_EXPONVAR_H_
#define SRC_RNG_EXPONVAR_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "Ran.h"

using namespace std;

/**
 * Structure for exponential variate.
 */
struct ExponVar : Ran {
	double beta;

	/** Constructor arguments are beta and a random sequence seed. */
	ExponVar(double beta): beta(beta) {}

	/** Returns an exponential deviate. */
	double var() {
		double u;
		do u = doub(); while (u == 0.);
		return -log(u)/beta;
	}
};

/**
 * Informal test of <ExponVar>. Prints readable test cases and results in the console.
 */
void test_ExponVar() {
	cout << "Testing Expondev." << endl << "-----------------" << endl;

	double beta = 1.5;
	ExponVar expon(beta);

	// Filling a vector <a> with <N> random entries.
	int N = 10000;
	double* a = new double[N];
	for(int i = 0; i < N; i++) a[i] = expon.var();

	// Initialization of quant array.
	int nQuant = 10;
	int* quant = new int[nQuant];
	for(int j = 0; j < nQuant; j++) quant[j] = 0;

	// Calculating the nQuant-quantiles of <a>.
	for(int j = 0; j < nQuant; j++){
		for(int i = 0; i < N; i++) {
			if (a[i] < -log(1 - (j+1)/(double)nQuant)/beta) quant[j] += 1;
		}
	}

	// Printing the results.
	cout << "Using the first " << N << " elements of the Expondev("<<beta<<") sequence, we have the following" << endl;
	cout << "estimated "<< nQuant<<"-quantiles (first line) versus the frequencies expected of exponential" << endl;
	cout << "random variables (second line). The seed is " << Ran::seed << " ." << endl;
	for(int j = 0; j < nQuant-1; j++) {
		cout << quant[j]/(double)N << "  ";
	}
	cout << endl;
	for(int j = 0; j < nQuant-1; j++) {
			cout << (j+1)/(double)nQuant << "  ";
		}
	cout << endl << endl;
}

#endif /* SRC_RNG_EXPONVAR_H_ */
