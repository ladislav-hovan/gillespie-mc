#include <iostream>
#include <fstream>

#include "input.h"
#include "gillespie.h"
#include "rungekutta.h"

using std::vector;
using std::pair;

template <typename T>
void printCountsToStream(T &Stream, countMap &Counts, bool bLabels=false)
{
	for (auto &Pair: Counts)
	{
		if (bLabels)
			Stream << Pair.first << " ";
		else
			Stream << Pair.second << " ";
	}
	Stream << "\n";
}

int main()
{
	// Obtain data about possible transitions and initial conditions
	Input sInput;
	vector<Reaction> vsReactions = loadInput("input.dat", sInput);
	if (vsReactions.empty())
	{
		std::cout << "No valid reactions were entered" << std::endl;
		std::getchar();
		exit(0);
	}

	// Map of initial counts
	countMap Counts;
	for (auto &Pair: sInput.initial)
		Counts[Pair.first] = Pair.second;

	// Keep track of allowed reactions (enough reactants)
	vector<bool> vbAllowed;
	for (auto &Reaction: vsReactions)
		vbAllowed.push_back(checkSufficient(Counts, Reaction));

	// Map of reactants to set of reactions (represented by position)
	rxnMap reactantsToRxns;
	for (int nPos = 0; nPos < vsReactions.size(); ++nPos)
	{
		for (auto &Pair: vsReactions[nPos].reactants)
		{
			// Create an empty set if the reactant isn't used yet
			reactantsToRxns.emplace(std::make_pair(Pair.first, std::set<unsigned int>()));
			// Add reaction to the set if it isn't there
			reactantsToRxns[Pair.first].insert(nPos);
		}
	}

	// Vectors of reaction rates and their cumulative sums
	vector<double> vdRates;
	for (int nPos = 0; nPos < vsReactions.size(); ++nPos)
		vdRates.push_back(getRate(vbAllowed, Counts, vsReactions, nPos));
	vector<double> vdCumSums = vdRates;
	if (!sInput.first)
	{
		for (int nPos = 1; nPos < vdCumSums.size(); ++nPos)
			vdCumSums[nPos] += vdCumSums[nPos - 1];
	}
	vector<double> &vdRatesOrSums = (sInput.first ? vdRates : vdCumSums);

	// Output file
	std::ofstream Output;
	if (sInput.outfile != "")
		Output.open(sInput.outfile);

	// Algorithm loop
	double dTime = 0.0;
	Output << "# Time ";
	printCountsToStream(Output, Counts, true);
	Output << dTime << " ";
	printCountsToStream(Output, Counts);
	std::cout << std::endl;
	for (int nStep = 0; nStep < sInput.steps; ++nStep)
	{
		auto Pair = iterateGillespie(vdRatesOrSums, sInput.first, sInput.binary);
		adjustEverything(Pair.first, vsReactions, Counts, vbAllowed, reactantsToRxns, vdRates,
			vdCumSums, sInput.first);
		dTime += Pair.second;
		if (sInput.outfile != "" && nStep % sInput.fprogress == sInput.fprogress - 1)
		{
			Output << dTime << " ";
			printCountsToStream(Output, Counts);
		}
		if (nStep % sInput.progress == sInput.progress - 1)
		{
			std::cout << "Step " << nStep + 1 << ", time " << dTime << ": ";
			printCountsToStream(std::cout, Counts);
		}
	}
	if (Output.is_open())
		Output.close();

	if (sInput.estimate)
	{
		std::cout << std::endl;
		countMapDouble Estimate = estimateEquilibrium(sInput.initial, reactantsToRxns, vsReactions,
			sInput.maxsteps, sInput.stepsize);
		std::cout << "Equilibrium estimates:" << std::endl;
		for (auto &Pair: Estimate)
			std::cout << Pair.first << ": " << std::round(Pair.second) << std::endl;
	}

	std::getchar();
	return 0;
}