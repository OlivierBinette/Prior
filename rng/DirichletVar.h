/*
 * DirichletVar.h
 *
 *  Created on: Dec 24, 2016
 *      Author: Binette
 */

#ifndef SRC_RNG_DIRICHLETVAR_H_
#define SRC_RNG_DIRICHLETVAR_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "GammaVar.h"

using namespace std;

/**
 * Structure for Dirichlet variates.
 *
 * USAGE
 * =====
 * >>> int N = 10; double* write_in = new double[N];	// Array in which to return the Dirichlet variates.
 * >>> double* alpha = new double[N]; 	// Positive parameters of the Dirichlet distribution.
 * >>> for (int i = 0; i < N; i++) alpha[i] =
 *
 */
struct DirichletVar {

	int size;
	double* alpha;
	vector<GammaVar> gammaVar;

	/**
	 * Constructor argument <aalpha> is the array of parameters for the distribution. Its entries
	 * must be strictly positive.
	 */
	DirichletVar(double* aalpha, int ssize) {
		size = ssize;
		alpha = new double[size];
		for (int i = 0; i < size; i++)  {
			if (aalpha[i] <= 0) throw invalid_argument("Non-positive parameter in DirichletVar.");
			alpha[i] = aalpha[i];
		}
		for (int i = 0; i < size; i++) gammaVar.push_back(GammaVar(alpha[i], 1));
	}

	void var(double write_in[]) {
		for (int i = 0; i < size; i++) write_in[i] = gammaVar[i].var();
		double s = sum(write_in, size);
		for (int i = 0; i < size; i++) write_in[i] /= s;
	}
};

#endif /* SRC_RNG_DIRICHLETVAR_H_ */
