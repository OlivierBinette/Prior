#ifndef SRC_RNG_NORMALVAR_H_
#define SRC_RNG_NORMALVAR_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "Ran.h"

using namespace std;

/**
 * Structure for normal variates (Numerical Recipes, 3rd ed., p.369).
 */
struct NormalVar : Ran {
	double mu, sig;

	/** Constructor arguments are mu, sigma and a random sequence seed. */
	NormalVar(double mu, double sig): mu(mu), sig(sig) {}

	NormalVar(const NormalVar& o) : Ran(o) {
		mu = o.mu; sig = o.sig;
	}

	/** Returns a normal deviate. */
	double var() {
		double u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = fabs(v) + 0.386595;
			q = x*x + y*(0.19600*y - 0.25472*x);
		} while (q > 0.27597 && (q > 0.27846 || v*v > -4.*log(u)*u*u));
		return mu + sig*v/u;
	}
};


#endif /* SRC_RNG_NORMALVAR_H_ */
