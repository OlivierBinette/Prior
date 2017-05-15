#ifndef SRC_RNG_RAN_H_
#define SRC_RNG_RAN_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

typedef unsigned long long Ullong;

using namespace std;


/**
 * Random number generator. This is an adaptation of Numerical Recipe's "suspenders-and-belts, full-body-armor,
 * never-any-doubt generator" (third edition, p.342).
 *
 * Constructor is called with an integer seed. The member functions <int64>, <doub>, and <int32>
 * return the next values in the random sequence, as a variable type indicated by their names. The
 * period of the generator is +/- 3.138 * 10^{57} .
 *
 * USAGE
 * =====
 * >>> Ran rnd;		// Instance of the base random number generator struct. It is automatically
 * 					// given a unique seed.
 * >>> rnd.doub()   // Returns a uniform variate in the half-open interval [0, 1).
 */
struct Ran {
	static Ullong seed; // Seed of the previous instance of Ran. It is incremented at each instanciation.

	Ullong u,v,w;

	/** Constructor.*/
	Ran(): v(4101842887655102017LL), w(1) {
		seed += 1;
		u = seed ^ v; int64();
		v = u; int64();
		w = v; int64();
	}

	Ran(const Ran& o) {
		u = o.u; v = o.v; w = o.w;
	}

	/** Returns 64-bit random integer. */
	inline Ullong int64() {
		u = u*2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}

	/** Returns random double-precision floating value in the range 0. to 1. */
	inline double doub() { return 5.42101086242752217E-20 * int64(); }

	/** Returns 32-bit random integer. */
	inline unsigned int int32() { return (unsigned int)int64(); }
};
Ullong Ran::seed = 0;

/**
 * Informal test of <Ran>. Prints readable test cases and results in the console.
 */
void test_Ran() {
	cout << "Testing Ran." << endl << "------------" << endl;
	Ran rnd;

	// Filling a vector <a> with <N> random entries.
	int N = 10000;
	double* a = new double[N];
	for(int i = 0; i < N; i++) a[i] = rnd.doub();

	// Printing the first numbers of Ran(42).
	cout << "Here are the first 50 numbers of the sequence Ran() with seed "<< Ran::seed <<" : " << endl;
	for(int k = 0; k < 50; k++) {
		cout << setprecision(5) << fixed << a[k] << "  ";
		if (k % 10 == 9) cout << endl;
	}
	cout << endl << endl;

	// Initialization of quant array.
	int nBin = 10;
	int* bin = new int[nBin];
	for (int j = 0; j < nBin; j++) bin[j] = 0;

	// Calculating the nQuant-quantiles of <a>.
	for(int j = 0; j < nBin; j++){
		for(int i = 0; i < N; i++) {
			if (a[i] < (j+1)/(double)nBin) bin[j] += 1;
		}
	}

	// Printing the results.
	cout << "Using the first " << N << " elements of the Ran() sequence, we have the following " << endl;
	cout << "uniform "<< nBin<<"-quantiles frequencies (first line) versus the frequencies expected of uniform" << endl;
	cout << "random variables (second line). The seed is " << Ran::seed << " ." << endl;
	for(int j = 0; j < nBin-1; j++) {
		cout << bin[j]/(double)N << "  ";
	}
	cout << endl;
	for(int j = 0; j < nBin-1; j++) {
			cout << (j+1)/(double)nBin << "  ";
		}
	cout << endl << endl;
}

#endif /* SRC_RNG_RAN_H_ */
