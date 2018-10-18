#include "gillespie.h"

inline bool isNearlyEqual(double dX, double dY)
{
	const double dEpsilon = 1e-8;
	return std::abs(dX - dY) <= dEpsilon * std::abs(dX);
}

double getRand()
{
	// PRNG, static so only created once
	static std::random_device Seed;
	static std::mt19937 Mersenne(Seed());
	static std::uniform_real_distribution<double> Rand(0.0, 1.0);

	return Rand(Mersenne);
}

bool checkSufficient(countMap &Counts, Reaction &sReaction)
{
	bool bSufficient = true;
	for (auto &Pair : sReaction.reactants)
	{
		if (Counts[Pair.first] < Pair.second)
			bSufficient = false;
	}

	return bSufficient;
}

void adjustEverything(int nPos, vector<Reaction> &vsReactions, countMap &Counts,
	vector<bool> &vbAllowed, rxnMap &reactantsToRxns, vector<double> &vdRates,
	vector<double> &vdCumSums, bool bFirst)
{
	std::set<unsigned int> affectedRxns;
	for (auto &Pair: vsReactions[nPos].reactants)
	{
		// Update amount
		Counts[Pair.first] -= Pair.second;
		// Add to the set of affected reactions
		for (auto &nRxn: reactantsToRxns[Pair.first])
			affectedRxns.insert(nRxn);
	}
	for (auto &Pair: vsReactions[nPos].products)
	{
		// Update amount
		Counts[Pair.first] += Pair.second;
		// Add to the set of affected reactions
		for (auto &nRxn: reactantsToRxns[Pair.first])
			affectedRxns.insert(nRxn);
	}
	for (auto &nPos: affectedRxns)
	{
		// Recheck sufficiency of reactants
		vbAllowed[nPos] = checkSufficient(Counts, vsReactions[nPos]);
		// Recalculate rates
		vdRates[nPos] = getRate(vbAllowed, Counts, vsReactions, nPos);
	}
	if (!bFirst)
	{
		for (int nPos = *(affectedRxns.begin()); nPos < vdCumSums.size(); ++nPos)
		{
			if (nPos == 0)
				vdCumSums[0] = vdRates[0];
			else
				vdCumSums[nPos] = vdCumSums[nPos - 1] + vdRates[nPos];
		}
	}
}

pair<int, double> iterateGillespieDirect(vector<double> &vdCumSums, bool bBinary)
{
	double dFinalSum = vdCumSums[vdCumSums.size() - 1];

	// Check in case no transitions possible
	if (isNearlyEqual(dFinalSum, 0.0))
	{
		std::cout << "No reactions are possible, quitting" << std::endl;
		std::getchar();
		exit(0);
	}

	int nWhich = 0;
	// Decide on transition
	if (vdCumSums.size() > 1)
	{
		// Else only one option, so no need to generate random number
		double dProduct = getRand() * dFinalSum;
		if (bBinary)
			// Binary search, logarithmic scaling
			nWhich = std::lower_bound(vdCumSums.begin(), vdCumSums.end(), dProduct) - vdCumSums.begin();
		else
			// Searching in order, linear scaling, efficient for small number of possibilities
			while (vdCumSums[nWhich] < dProduct) ++nWhich;
	}

	// Calculate tau (time increment)
	double dTau = -(1 / dFinalSum) * log(getRand());

	return pair<int, double>(nWhich, dTau);
}

pair<int, double> iterateGillespieFirst(vector<double> &vdRates)
{
	// Calculate tau (time increment) for each possible target, keep smallest
	int nSmallest = -1;
	double dSmallest = std::numeric_limits<double>::infinity();
	for (unsigned int nPos = 0; nPos < vdRates.size(); ++nPos)
	{
		double dRate = vdRates[nPos];
		// A rate of zero can't enter into the next equation, so just ignore
		if (isNearlyEqual(dRate, 0.0))
			continue;
		double dTau = (1 / dRate) * log(1 / getRand());
		if (dTau < dSmallest)
		{
			nSmallest = nPos;
			dSmallest = dTau;
		}
	}
	if (nSmallest == -1)
	{
		std::cout << "No reactions are possible, quitting" << std::endl;
		std::getchar();
		exit(0);
	}

	// Go with the smallest tau
	return pair<int, double>(nSmallest, dSmallest);
}

pair<int, double> iterateGillespie(vector<double> &vdRatesOrSums, bool bFirst, bool bBinary)
{
	if (bFirst)
		return iterateGillespieFirst(vdRatesOrSums);
	else
		return iterateGillespieDirect(vdRatesOrSums, bBinary);
}