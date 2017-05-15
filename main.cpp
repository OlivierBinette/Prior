/*
 * main.cpp
 *
 *  Created on: May 3, 2017
 *      Author: Binette
 */

#include "prior.hpp"
#include "sampler.hpp"
#include <math.h>
#include <iostream>
#include "utils.hpp"
#include "parser.hpp"

using namespace std;

int main(int argc, char* argv[]) {

	/*
	Poussin basis;
	Poussin::UNIFORM_DISTRIB unif;

	Prior<double> p(basis, unif, 1, 5);

	vector<double> data(50, M_PI);
	data.push_back(0);
	data.push_back(1);
	data.push_back(1);
	data.push_back(1);
	data.push_back(1);
	data.push_back(1);
	data.push_back(1);
	data.push_back(1);



	p.update(data);


	Sampler<double> s(p);

	s.sample(1000000, 20);

	//vector<double>& mean = o.mean();

	//cout << endl << endl;
	//print_array(mean, mean.size(), 5);
	//cout << endl;

	//print_array(s.alpha(4), s.alpha(4).size(), 5);

	//cout << endl << endl;

	//print_array(s.alpha(5), s.alpha(5).size(), 5);

	//cout << endl<< endl;

	//print_array(s.alpha(6), s.alpha(6).size(), 5);

	vector<double> Y(11, 1./11);
	vector<double> util;

	//cout << endl << endl << endl << endl;


	// Parsing JSON input

	/*
	using json = nlohmann::json;
	string input = "";
	for (int i = 1; i < argc; i++)  input += argv[i];
	json j = json::parse(input);
	vector<double> data;
	for (int i =0; i < j["data"].size(); i++) data.push_back(j["data"][i]);
	double spread = j["spread"];
	double mu = j["mu"];
	int n_0 = j["n_0"];

	Poussin basis;
	Poussin::UNIFORM_DISTRIB unif;
	Prior<double> p(basis, unif, spread, mu);
	Sampler<double> s(p);
	p.update(data);
	Output o = s.sample(j["n_samples"], n_0);

	cout << o.to_json();
	*/

	string input = "";
	for (int i = 1; i < argc; i++)  input += argv[i];

	Parser params; params.run(input);
	Poussin basis;
	Poussin::UNIFORM_DISTRIB unif;

	Prior<double> p(basis, unif, params.concentration, params.smoothness);
	p.update(params.data);
	Sampler<double> s(p);
	s.max_dim = params.max_dim;
	s.lowest_dim = params.lowest_dim;

	s.sample(params.n_samples, params.n_0);

	return 0;
}

