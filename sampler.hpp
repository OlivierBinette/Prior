
#ifndef BINETTE_SAMPLER_HPP_
#define BINETTE_SAMPLER_HPP_

#include <iostream>
#include <vector>
#include "prior.hpp"

using namespace std;

struct Output {
	vector<vector<double> > samples;
	vector<double> _mean;

	Output() {}

	void append(vector<double> sample) {
		samples.push_back(sample);
	}

	// TODO Temporary definition
	vector<double>& mean() {
		double s;
		for (int j = 0; j < samples[0].size(); j++){
			s = 0;
			for (int i = 0; i < samples.size(); i++)
				s += samples[i][j];
			_mean.push_back(s/samples.size());
		}
		return _mean;
	}; //TODO

	string to_json() {
		return array_to_string(samples);
	}

	vector<vector<double> > mean_trajectory() {}; //TODO
};

template <typename X>
struct Sampler {
	Prior<X>& P;
	Ran rng;

	// Simulation parameters.
	double p_change_dim = 0.1;
	int lowest_dim = 1;
	double max_dim = 200;

	// Precalculated values
	vector<vector<double> > _alpha;
	vector<double> _lambda; // Sum of the _alpha rows.
	vector<vector<double> > _prop_alpha;
	vector<vector<vector<double> > > _values;
	vector<double> _poisson_log_pmf;

	Sampler(Prior<X>& P): P(P) {
		_alpha.push_back(vector<double>());
		_lambda.push_back(NAN);
		_prop_alpha.push_back(vector<double>());
		_values.push_back(vector<vector<double> >());
		_poisson_log_pmf.push_back(-P.mu);
	}

	Output sample(int n_times, int n_0) {
		Output out;

		double accepted = 0;
		double total = 0;
		bool moved = false;

		// Initialization
		int m_n_0 = P.basis.m(n_0);
		vector<double> Y(m_n_0, 1./m_n_0);
		vector<double> proposal(m_n_0, 0);
		vector<unsigned int> comps(2, 0);
		vector<double> util_vect; //
		int n = n_0; // Current degree/dimension.
		int m_N = P.basis.m(n);

		double A, u_1, u_2;
		do {
			// Propose new dimension.
			int prop_n = q_J(n);

			if (n == prop_n) {
				// Gibbs update of k components of Y.

				// New proposal
				for (int i = 0; i < Y.size(); i++) proposal[i] = Y[i];
				chooseComps(comps, 2, P.basis.m(n));
				for (int i = 0; i < 2; i++)
					proposal[comps[i]] = gamma_variate(prop_alpha(n)[comps[i]], rng);

				// Acceptance probability.
				A = log_posterior_pdf(proposal, util_vect, n, false) + q_log_pdf(Y, comps, 2, n)
						- log_posterior_pdf(Y, util_vect, n, false) - q_log_pdf(proposal, comps, 2, n);

				//cout << "Gibbs step. A = " << exp(A) << endl;

				if (log(rng.doub()) < A){
					//cout << "Updated. n = " << prop_n << endl;
					accepted += 1;
					moved = true;
					for (int i = 0; i < 2; i++)
						Y[comps[i]] = proposal[comps[i]];
				}
			}
			else if (n < prop_n) { // prop_n = n+1
				// Split move.

				m_N = P.basis.m(n);
				u_1 = gamma_variate(prop_alpha(prop_n)[m_N], rng);
				u_2 = gamma_variate(prop_alpha(prop_n)[m_N + 1], rng);

				proposal.resize(Y.size() + 2);
				for (int i = 0; i < Y.size(); i++) proposal[i] = Y[i];
				proposal[m_N] = u_1;
				proposal[m_N + 1] = u_2;


				A = log_posterior_pdf(proposal, util_vect, prop_n, true) - log_posterior_pdf(Y, util_vect, n, true);
				A += log_J(n, prop_n) - log_J(prop_n, n);
				A -= gamma_log_pdf(u_1, prop_alpha(prop_n)[m_N]);
				A -= gamma_log_pdf(u_2, prop_alpha(prop_n)[m_N + 1]);

				//cout << setprecision(6) << "Splitting. A = " << exp(A) << ". u_1 = " << u_1 << ". u_2 = " << u_2 <<"." << endl;

				if (log(rng.doub()) < A) {
					//cout << "Splitted. n = " << prop_n << endl;
					//cout << "===========================================================================================================" << endl;
					accepted += 1;
					moved = true;
					Y.resize(Y.size() + 2);
					Y[m_N] = u_1;
					Y[m_N+1] = u_2;
					n = prop_n;
				}
			}
			else if (n > prop_n) {
				// Merge move.
				m_N = P.basis.m(n);
				proposal.resize(Y.size()-2);
				for(int i = 0; i < Y.size()-2; i++) proposal[i] = Y[i];
				u_1 = Y[m_N-2];
				u_2 = Y[m_N-1];

				A = log_posterior_pdf(proposal, util_vect, prop_n, true) - log_posterior_pdf(Y, util_vect, n, true);
				A += log_J(n, prop_n) - log_J(prop_n, n);

				/*cout << "log_post(prop): " << log_posterior_pdf(proposal, util_vect, prop_n, true) << endl;
				cout << "log_post(Y): " << log_posterior_pdf(Y, util_vect, prop_n, true) << endl;
				cout << "log_J(n, prop_n) :" << log_J(n, prop_n) << endl;
				cout << "log_J(prop_n, n) :" << log_J(prop_n, n) << endl;
				*/
				A += gamma_log_pdf(u_1, prop_alpha(n)[m_N - 2]);
				A += gamma_log_pdf(u_2, prop_alpha(n)[m_N - 1]);

				//cout << "gamma_log_pdf (u_1) : " << gamma_log_pdf(u_1, prop_alpha(prop_n)[m_N-2]) << endl;
				//cout << "gamma_log_pdf (u_2) : " << gamma_log_pdf(u_1, prop_alpha(prop_n)[m_N-1]) << endl;

				//cout << setprecision(6) << "Merging. A = " << A << ". u_1 = " << u_1 << ". u_2 = " << u_2 <<"." << endl;

				if (log(rng.doub()) < A) {
					//cout << "Merged. n = "<< prop_n << endl;
					//cout << "===========================================================================================================" << endl;
					accepted += 1;
					moved = true;
					Y.resize(Y.size()-2);
					n = prop_n;
				}
			}
			total += 1;
			if (moved) {
				phiTransform(util_vect, Y);
				util_vect.resize(Y.size());
				out.append(util_vect);
			}
			moved = false;
			//print_array(util_vect, P.basis.m(n), 3);
			//cout << endl;

			n_times -= 1;
		} while (n_times > 0);

		//cout << endl << "Acceptance ratio : " << 100. * accepted/total << " %" << endl << endl;

		//vector<double> mean = out.mean();
		//print_array(mean, mean.size(), 5);
		cout << out.to_json() << endl;
	}

	inline double log_posterior_pdf(vector<double>& Y, vector<double>& c_theta, int n, bool consts) {
		phiTransform(c_theta, Y);
		double p = 0, m_N = P.basis.m(n);

		// Log likelihood.
		for (int i = 0; i < values(n).size(); i++) {
			p += log(dot(c_theta, values(n)[i], m_N));
		}

		// Log prior.
		for (int i = 0; i < m_N; i++) {
			if (consts){
				p += (alpha(n)[i]-1) * log(c_theta[i]) - lgamma(alpha(n)[i]); // Dirichlet kernel.
			}
			else {
				p += (alpha(n)[i]-1) * log(c_theta[i]);
			}
			// The gamma(sum(a_i)) constant cancels out in the Gamma(theta, lambda(n) density.
		}
		p += (lambda(n)-1)*log(c_theta[m_N]) - c_theta[m_N]; // Gamma(theta, lambda(n)) density

		p += poisson_log_pmf(n);


		return p + (1-m_N)*log(sum(Y, m_N));
	}

	inline double q_log_pdf(const vector<double>& Y, const vector<unsigned int>& comps, int k, int n) {
		double s = 0.0;
		for (int i = 0; i < k; i++) {
			s += (prop_alpha(n)[comps[i]] - 1) * log(Y[comps[i]]) - Y[comps[i]];
		}
		return s;
	}

	/**
	 * Returns the nth degree prior parameters for the density coefficients.
	 */
	vector<double>& alpha(int n) {
		if (n < _alpha.size()) return _alpha[n];
		alpha(n-1);
		_alpha.push_back(P.basis.proj(P.center, n, P.spread));
		return _alpha[n];
	}

	/**
	 * Returns the sum of the nth degree prior parameters for the density coefficients.
	 */
	double lambda(int n) {
		if (n < _lambda.size()) return _lambda[n];
		lambda(n-1);
		_lambda.push_back(sum(alpha(n), alpha(n).size()));
		return _lambda[n];
	}

	/*
	 * Returns the nth degree proposal parameters for the density coefficients.
	 */
	vector<double>& prop_alpha(int n) {
		if (n < _prop_alpha.size()) return _prop_alpha[n];
		prop_alpha(n-1);

		vector<int> bins = P.basis.bin(P.data, n);
		vector<double> prop_a = vector<double>(P.basis.m(n), 0);
		vector<double>& alph = alpha(n);
		for (int i = 0; i < P.basis.m(n); i++) prop_a[i] = alph[i] * (1. + bins[i]);

		_prop_alpha.push_back(prop_a);
		return _prop_alpha[n];
	}

	/**
	 * Returns the values of the nth degree densities at the data points.
	 *
	 * values[i] is the vector of the densities evaluated at the data point x_i.
	 * values[i][j] is the jth density evaluated at x_i.
	 *
	 */
	vector<vector<double> >& values(int n) {
		if (n < _values.size()) return _values[n];
		values(n-1);

		vector<vector<double> > vals;
		for (int i = 0; i < P.data.size(); i++) {
			vals.push_back(vector<double>());
			for (int j = 0; j < P.basis.m(n); j++) {
				vals[i].push_back(P.basis.density(P.data[i], j, n));
			}
		}

		_values.push_back(vals);
		return _values[n];
	}

	double log_J(int prop_n, int n) {
		if (prop_n == n) return log(1-p_change_dim);
		if (prop_n < n) return log(p_change_dim *((n-lowest_dim) / (max_dim - lowest_dim)));
		return log(p_change_dim * (1 - (n-lowest_dim) / (max_dim - lowest_dim)));
	};

	int q_J(int n) {
		double u = rng.doub();
		if (u <= 1-p_change_dim) return n;
		u -= 1-p_change_dim;
		if (u < p_change_dim*(1 - (n-lowest_dim)/(max_dim - lowest_dim)) && n < max_dim) return n+1;
		return n-1;
	};

	inline void phiTransform(vector<double>& write_in, const vector<double>& Y) {
		write_in.resize(Y.size() + 1);
		double s = sum(Y, Y.size());
		for (int i = 0; i < Y.size(); i++) write_in[i] = Y[i] / s;
		write_in[Y.size()] = s;
	};

	inline void chooseComps(vector<unsigned int>& write_in, int k, int n) {
		double prop; bool unique;
		write_in[0] = rng.int32() % n;
		for (int i = 1; i < k; i++) {
			do {
				unique = true;
				prop = rng.int32() % n;
				for (int j = 0; j < i; j++) {
					if (write_in[j] == prop) unique = false;
				}
			} while (!unique);
			write_in[i] = prop;
		}
	}

	double poisson_log_pmf(int n) {
		if (n < _poisson_log_pmf.size()) return _poisson_log_pmf[n];
		poisson_log_pmf(n-1);

		double log_n_factorial = 0.0;
		double c = 0.0, y, t;
		for (int i = 1; i <= n; i++) {
			y = log(i) - c;
			t = log_n_factorial + y;
			c = (t - log_n_factorial) - y;
			log_n_factorial = t;
		}
		_poisson_log_pmf.push_back(-P.mu + n*log(P.mu) - log_n_factorial);
		return _poisson_log_pmf[n];
		//return -n;
	}
};


#endif
