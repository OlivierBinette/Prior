/*
 * utils.hpp
 *
 *  Created on: May 3, 2017
 *      Author: Binette
 */

#ifndef BINETTE_UTILS_HPP_
#define BINETTE_UTILS_HPP_

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <sstream>
#include "rng/Ran.h"

using namespace std;

/**
 * Returns the sum of the values of the array <a>, using Kahan's summation algorithm.
 */
double sum(const vector<double>& a, int size) {
	double s = 0.0, c = 0.0, y, t;
	for (int i = 0; i < size; i++) {
		y = a[i] - c;
		t = s + y;
		c = (t - s) - y;
		s = t;
	}
	return s;
}

double sum(double a[], int size) {
	double s = 0.0, c = 0.0, y, t;
	for (int i = 0; i < size; i++) {
		y = a[i] - c;
		t = s + y;
		c = (t - s) - y;
		s = t;
	}
	return s;
}

void print_array(vector<double>& a, int size, int n_decimals) {
	cout << "[ ";
	for (int i = 0; i < size; i++) {
		cout << setprecision(n_decimals) << fixed << a[i] << " ";
	}
	cout << "]";
}

string array_to_string(vector<vector<double> >& a) {
	ostringstream s;

	s << "[";
	for (int i = 0; i < a.size(); i++){
		s << "[";
		for (int j = 0; j < a[i].size(); j++) {
			s << a[i][j];
			if (j < a[i].size() - 1) s << ", ";
		}
		s << "]";
		if (i < a.size() - 1) s << ", ";
	}
	s << "]";
	return s.str();
}

/**
 * Returns the dot product of two vectors.
 */
double dot(vector<double>& a, vector<double>& b, int size) {
	double s = 0.0;
	for (int i = 0; i < size; i++) s += a[i]*b[i];
	return s;
}

/***
 * Statistical utilities
 */

inline double std_normal_variate(Ran& rng) {
	double u,v,x,y,q;
	do {
		u = rng.doub();
		v = 1.7156*(rng.doub()-0.5);
		x = u - 0.449871;
		y = fabs(v) + 0.386595;
		q = x*x + y*(0.19600*y - 0.25472*x);
	} while (q > 0.27597 && (q > 0.27846 || v*v > -4.*log(u)*u*u));
	return v/u;
}

/**
 * Gamma(alpha, 1) variate.
 */
double gamma_variate(double alpha, Ran& rng) {
	double oalpha = alpha;
	if (alpha < 1.) alpha += 1.;
	double a1 = alpha - 1./3.;
	double a2 = 1./sqrt(9.*a1);


	double u,v,x;
	do {
		do {
			x = std_normal_variate(rng);
			v = 1. + a2*x;
		} while (v <= 0.);
		v = v*v*v;
		u = rng.doub();
	} while(u > 1. - 0.331*x*x*x*x && log(u) > 0.5*x*x + a1*(1. - v + log(v)));
	if (alpha == oalpha) return a1*v;
	else {
		do u = rng.doub(); while (u == 0.);
		return pow(u, 1./oalpha)*a1*v;
	}
}

double gamma_log_pdf(double u, double alpha) {
	return (alpha-1)*log(u) - lgamma(alpha) - u;
};

#endif /* SRC_UTILS_HPP_ */
