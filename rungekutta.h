#pragma once

#include "reaction.h"
#include "gillespie.h"

using countMapDouble = std::unordered_map<unsigned int, double>;

double calculateDerivative(int nID, double dCount, rxnMap &reactantsToRxns, rxnMap &productsToRxns,
	vector<Reaction> &vsReactions, countMapDouble &Counts);
countMapDouble iterateRungeKutta(countMapDouble &Counts, rxnMap &reactantsToRxns,
	vector<Reaction> &vsReactions, double dTimeStep = 0.1);
countMapDouble estimateEquilibrium(vector<pair<unsigned int, unsigned int>> &vpnInitial,
	rxnMap &reactantsToRxns, vector<Reaction> &vsReactions,
	unsigned int nMaxSteps = 1000, double dTimeStep = 0.1);