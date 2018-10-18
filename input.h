#pragma once

#include <fstream>
#include <sstream>
#include <cctype>
#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

#include "reaction.h"

using std::vector;
using std::pair;
using std::string;

struct Input
{
	bool first = false;
	bool binary = false;
	vector<pair<unsigned int, unsigned int>> initial;
	int steps = 100;
	int progress = 100;
	int fprogress = 1;
	string outfile = "";
	// The following concerns Runge-Kutta integration
	bool estimate = false;
	int maxsteps = 1000;
	double stepsize = 0.1;
};

void registerReaction(Reaction &sRxn, vector<string> &vsWords, Input &sInput);
Reaction analyseLine(string strLine, Input &sInput);
vector<Reaction> loadInput(string strFilename, Input &sInput);