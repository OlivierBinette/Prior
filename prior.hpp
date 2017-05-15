/*
 * Prior. A statistical package for semiparametric mixture models.
 *
 *
 *
 *
 */

#ifndef BINETTE_PRIOR_HPP
#define BINETTE_PRIOR_HPP

#include <vector>
#include <math.h>
#include <iostream>
#include <functional>
#include "utils.hpp"

using namespace std;
using namespace std::placeholders;

/**
 * Probability distribution on the sample space X.
 */
template <typename X>
struct Distribution {

	/**
	 * Probability density.
	 */
	virtual double operator ()(X x) const {};
};

/**
 * The density basis \{ \phi_{i,n} \}_{i=1}^{m_n} with its associated partition \{ R_{i,n} \}_{i=1}^{m_n}
 * of the sample space. Each R_{i,n} is referred to as the region associated to \phi_{i,n}.
 */
template <typename X>
struct Basis {
	virtual ~Basis() {};

	/**
	 * The number of regions densities for a given degree n.
	 */
	virtual int m(int n) {};

	/**
	 * The i'th density of degree n, evaluated at x.
	 */
	virtual double density(X x, int i, int n) {};

	/**
	 * Projects the distribution d on the m(n) regions associated with the family of densities.
	 */
	virtual vector<double> proj(const Distribution<X>& d, int n, double mult_factor = 1.) {};

	/**
	 * Counts the number of data points falling in each of the m(n) regions associated with the
	 * family of densities.
	 */
	virtual vector<int> bin(const vector<X>& data, int n) {};
};

template <typename X>
struct Prior {
	Basis<X>& basis;
	Distribution<X>& center;
	double spread;
	double mu;

	vector<X> data;

	Prior (Basis<X>& basis, Distribution<X>& center, double spread, double mu):
		basis(basis), center(center), spread(spread), mu(mu) {}

	void update(const vector<X>& d) {
		data = d;
	}
};

/**
 * De la Vallee Poussin basis on the circle.
 */
class Poussin: public Basis<double> {
public:

	vector<double> normalization_constants;

	Poussin() {
		normalization_constants.push_back(1./(2*M_PI));
		c(200); // Precalculate the first 200 normalization constants.
	}

	int m(int n) {
		return 2*n + 1;
	}

	double density(double u, int j, int n) {
		return (2*n+1) * c(n) * pow((1 + cos(u - 2 * M_PI * j / (2*n+1.)))/2., n);
	}

	vector<double> proj(const Distribution<double>& d, int n, double mult_factor = 1.) {
		vector<double> values;
		for (int i = 0; i < m(n); i++) {
			values.push_back(mult_factor * integrate(d, M_PI*(2.*i-1.)/m(n), M_PI*(2.*i+1.)/m(n), 50));
		}

		return values;
	}

	vector<int> bin(const vector<double>& data, int n) {
		int n_bins = m(n);
		vector<int> bins(m(n), 0);

		for (int i=0; i < data.size(); i++) {
			if (data[i] >= M_PI*(2*n_bins - 1) / n_bins || data[i] < M_PI/n_bins) bins[0] += 1;
			else for (int j = 1; j < n_bins; j++) {
				if (data[i] >= M_PI*(2*j -1) / n_bins && data[i] < M_PI*(2*j+1)/n_bins) bins[j] += 1;
			}
		}

		return bins;
	}

	struct UNIFORM_DISTRIB: Distribution<double> {
		double operator()(double u) const {
			return 1./(2*M_PI);
		}
	};

	/**
	 * Returns the nth normalization constant. Values up to the largest n with which c was called
	 * are stored.
	 */
	inline double c(int n) {
		if (n < normalization_constants.size()) return normalization_constants[n];
		else {
			double s = c(n-1) * 2 * n /(2.*n+1.);
			normalization_constants.push_back(s);
			return s;
		}
	}

	double integrate(const Distribution<double>& d, double a, double b, int n) {
		double s = 0., c = 0., y, t;
		for (int i = 0; i < n; i++) {
			y = d(a + i*(b-a)/n) - c;
			t = s + y;
			c = (t - s) - y;
			s = t;
		}
		return (b-a) * s / n;
	}
};


#endif
