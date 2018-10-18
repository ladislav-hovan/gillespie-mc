#include "rungekutta.h"

double calculateDerivative(int nID, double dCount, rxnMap &reactantsToRxns, rxnMap &productsToRxns,
	vector<Reaction> &vsReactions, countMapDouble &Counts)
{
	static vector<bool> vbAllAllowed(vsReactions.size(), true);
	countMapDouble CountsCopy = Counts;
	CountsCopy[nID] = dCount;

	double dDerivative = 0.0;
	for (auto &nRxn: reactantsToRxns[nID])
	{
		int nCoeff = 0;
		for (auto &Pair: vsReactions[nRxn].reactants)
		{
			if (Pair.first == nID)
			{
				nCoeff = Pair.second;
				break;
			}
		}

		dDerivative -= (getRate(vbAllAllowed, CountsCopy, vsReactions, nRxn) * nCoeff);
	}
	for (auto &nRxn : productsToRxns[nID])
	{
		int nCoeff = 0;
		for (auto &Pair : vsReactions[nRxn].products)
		{
			if (Pair.first == nID)
			{
				nCoeff = Pair.second;
				break;
			}
		}

		dDerivative += (getRate(vbAllAllowed, CountsCopy, vsReactions, nRxn) * nCoeff);
	}

	return dDerivative;
}

countMapDouble iterateRungeKutta(countMapDouble &Counts, rxnMap &reactantsToRxns, 
	vector<Reaction> &vsReactions, double dTimeStep)
{
	// Map of products to set of reactions (represented by position)
	rxnMap productsToRxns;
	for (int nPos = 0; nPos < vsReactions.size(); ++nPos)
	{
		for (auto &Pair : vsReactions[nPos].products)
		{
			// Create an empty set if the product isn't used yet
			productsToRxns.emplace(std::make_pair(Pair.first, std::set<unsigned int>()));
			// Add reaction to the set if it isn't there
			productsToRxns[Pair.first].insert(nPos);
		}
	}

	vector<double> vdDeltas;
	for (auto &Pair: Counts)
	{
		vector<double> vdIntermediates(4);
		vdIntermediates[0] = calculateDerivative(Pair.first, Pair.second, 
			reactantsToRxns, productsToRxns, vsReactions, Counts);
		vdIntermediates[1] = calculateDerivative(Pair.first, Pair.second + (vdIntermediates[0] *
			dTimeStep / 2), reactantsToRxns, productsToRxns, vsReactions, Counts);
		vdIntermediates[2] = calculateDerivative(Pair.first, Pair.second + (vdIntermediates[1] *
			dTimeStep / 2), reactantsToRxns, productsToRxns, vsReactions, Counts);
		vdIntermediates[3] = calculateDerivative(Pair.first, Pair.second + (vdIntermediates[2] *
			dTimeStep), reactantsToRxns, productsToRxns, vsReactions, Counts);
		vdDeltas.push_back(dTimeStep * (vdIntermediates[0] + 2 * vdIntermediates[1] + 2 *
			vdIntermediates[2] + vdIntermediates[3]) / 6);
	}

	int nPos = 0;
	countMapDouble NewCounts;
	for (auto &Pair: Counts)
	{
		NewCounts[Pair.first] = Pair.second + vdDeltas[nPos++];
		if (NewCounts[Pair.first] < 0.0)
			NewCounts[Pair.first] = 0.0;
	}

	return NewCounts;
}

countMapDouble estimateEquilibrium(vector<pair<unsigned int, unsigned int>> &vpnInitial,
	rxnMap &reactantsToRxns, vector<Reaction> &vsReactions,
	unsigned int nMaxSteps, double dTimeStep)
{
	// Convert integer counts into doubles to facilitate integration
	countMapDouble Counts;
	for (auto &Pair: vpnInitial)
		Counts[Pair.first] = static_cast<double>(Pair.second);

	// Iterate until no change or max steps reached
	for (int nStep = 0; nStep < nMaxSteps; ++nStep)
	{
		countMapDouble NewCounts = iterateRungeKutta(Counts, reactantsToRxns, vsReactions, dTimeStep);
		bool bSame = true;
		for (auto &Pair: Counts)
		{
			if (!isNearlyEqual(Pair.second, NewCounts[Pair.first]))
				bSame = false;
		}

		if (bSame)
		{
			std::cout << "Runge-Kutta integration has converged" << std::endl;
			return NewCounts;  // Or Counts, it's the same thing
		}
		else
			Counts = NewCounts;
	}

	// Integration hasn't converged
	std::cout << "Runge-Kutta integration has not converged" << std::endl;
	return Counts;
}