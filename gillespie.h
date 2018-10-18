#pragma once

#include <unordered_map>
#include <set>
#include <vector>
#include <utility>
#include <random>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "reaction.h"

using countMap = std::unordered_map<unsigned int, unsigned int>;
using rxnMap = std::unordered_map<unsigned int, std::set<unsigned int>>;
using std::vector;
using std::pair;

inline bool isNearlyEqual(double dX, double dY);
double getRand();
bool checkSufficient(countMap &Counts, Reaction &sReaction);
void adjustEverything(int nPos, vector<Reaction> &vsReactions, countMap &Counts,
	vector<bool> &vbAllowed, rxnMap &reactantsToRxns, vector<double> &vdRates,
	vector<double> &vdCumSums, bool bFirst);
pair<int, double> iterateGillespieDirect(vector<double> &vdRates, bool bBinary);
pair<int, double> iterateGillespieFirst(vector<double> &vdRates);
pair<int, double> iterateGillespie(vector<double> &vdRates, bool bFirst = false,
	bool bBinary = false);

// This allows to use both integer and double values for counts
template <typename T>
double getRate(vector<bool> &vbAllowed, std::unordered_map<unsigned int, T> &Counts, 
	vector<Reaction> &vsReactions, int nPos)
{
	if (!vbAllowed[nPos])
		return 0.0;

	Reaction &sReaction = vsReactions[nPos];
	switch (sReaction.type)
	{
	case RxnType::Influx:
		return sReaction.c;
	case RxnType::Mono:
		return sReaction.c * Counts[sReaction.reactants[0].first];
	case RxnType::AB:
		return sReaction.c * Counts[sReaction.reactants[0].first] *
			Counts[sReaction.reactants[1].first];
	case RxnType::A2:
		return sReaction.c * Counts[sReaction.reactants[0].first] *
			(Counts[sReaction.reactants[0].first] - 1) / 2;
	case RxnType::ABC:
		return sReaction.c * Counts[sReaction.reactants[0].first] *
			Counts[sReaction.reactants[1].first] *
			Counts[sReaction.reactants[2].first];
	case RxnType::A3:
		return sReaction.c * Counts[sReaction.reactants[0].first] *
			(Counts[sReaction.reactants[0].first] - 1) *
			(Counts[sReaction.reactants[0].first] - 2) / 6;
	case RxnType::A2B:
		// It is guaranteed that the first one is A
		return sReaction.c * Counts[sReaction.reactants[0].first] *
			(Counts[sReaction.reactants[0].first] - 1) *
			Counts[sReaction.reactants[1].first] / 2;
	default:
		std::cout << "Reaction of unknown type chosen, terminating" << std::endl;
		std::getchar();
		exit(0);
	}
}