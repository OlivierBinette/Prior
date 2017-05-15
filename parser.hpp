/*
 * parser.hpp
 *
 *  Created on: May 9, 2017
 *      Author: Binette
 */

#include <regex>
#include <stdlib.h>
#include "utils.hpp"
#include <iostream>

using namespace std;

struct Parser {
	double lowest_dim = 1;
	int max_dim = 100;
	double p_change_dim = 0.3;
	double concentration = 1.;
	double smoothness = 3.;
	int n_samples = 10000;
	int n_0 = 5;
	vector<double> data;

	void run(string s) {

		string float_capture = "([0-9]+(\\.[0-9]+)?).*";
		string array_capture = "\\[(?:([0-9]+(\\.[0-9]+)?)\\,)*?([0-9]+(?:\\.[0-9]+)?)\\]";

		regex rlowest_dim (".*(?:lowest_dim:)"+float_capture);
		regex rmax_dim (".*(?:max_dim:)"+float_capture);
		regex rp_change_dim (".*(?:p_change_dim:)"+float_capture);
		regex rconcentration (".*(?:concentration:)"+float_capture);
		regex rsmoothness (".*(?:smoothness:)"+float_capture);
		regex rn_samples (".*(?:n_samples:)"+float_capture);
		regex rn_0 (".*(?:n_0:)"+float_capture);
		regex rdata (".*(?:data:)("+array_capture+").*");


		smatch sm;

		regex_match(s, sm, rlowest_dim);
		if (sm.size() > 0) lowest_dim = atof(sm[1].str().c_str()); //TODO: validate input
		regex_match(s, sm, rmax_dim);
		if (sm.size() > 0) max_dim = atof(sm[1].str().c_str()); //TODO: validate input
		regex_match(s, sm, rp_change_dim);
		if (sm.size() > 0) p_change_dim = atof(sm[1].str().c_str()); //TODO: validate input
		regex_match(s, sm, rconcentration);
		if (sm.size() > 0) concentration = atof(sm[1].str().c_str()); //TODO: validate input
		regex_match(s, sm, rsmoothness);
		if (sm.size() > 0) smoothness = atof(sm[1].str().c_str()); //TODO: validate input
		regex_match(s, sm, rn_samples);
		if (sm.size() > 0) n_samples = atof(sm[1].str().c_str()); //TODO: validate input
		regex_match(s, sm, rn_0);
		if (sm.size() > 0) n_0 = atof(sm[1].str().c_str()); //TODO: validate input

		regex_match(s, sm, rdata);
		string arr = sm[1].str();
		regex part("([0-9]+(\\.[0-9]+)?)\\,?");
		string n;
		while(regex_search(arr, sm, part))
		{
			n = sm.str();
			if (n[n.size()-1] == ',') n = n.substr(0, n.size()-1);
			data.push_back(atof(n.c_str()));
			arr = sm.suffix();
		}
	}

};


