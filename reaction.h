#pragma once

#include <vector>
#include <utility>

using std::vector;
using std::pair;

enum class RxnType
{
	Influx,  // No reactants
	Mono,	 // A -> products
	AB,		 // A + B -> products
	A2,		 // 2A -> products
	ABC,	 // A + B + C -> products
	A2B,	 // 2A + B -> products
	A3		 // 3A -> products
};

struct Reaction
{
	// Contains IDs and coefficients
	vector<pair<unsigned int, signed int>> reactants;
	vector<pair<unsigned int, signed int>> products;
	double c = -1.0;
	RxnType type;
};