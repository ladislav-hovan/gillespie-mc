#include "input.h"

void registerReaction(Reaction &sRxn, vector<string> &vsWords, Input &sInput)
{
	unsigned int nOrder = 0;
	bool bLeft = true;

	for (auto &strWord: vsWords)
	{
		if (strWord == "-")
		{
			if (bLeft)
				bLeft = false;
			else
				std::cout << "Ignoring additional \"-\" in reaction definition" << std::endl;
		}
		else if (strWord[0] == 'c')
		{
			string strVal = strWord.substr(2, strWord.size() - 2);
			sRxn.c = std::atof(strVal.c_str());
		}
		else if (std::isdigit(strWord[0]))
		{
			signed int nCoeff;
			unsigned int nID;

			int nPos = strWord.find('*');
			if (nPos < strWord.size())
			{
				nCoeff = abs(std::atoi(strWord.substr(0, nPos).c_str()));
				nID = std::atoi(strWord.substr(nPos + 1, strWord.size() - nPos - 1).c_str());
			}
			else
			{
				nCoeff = 1;
				nID = std::atoi(strWord.c_str());
			}

			vector<pair<unsigned int, signed int>>* side;
			if (bLeft)
			{
				nOrder += nCoeff;
				side = &(sRxn.reactants);
			}
			else
				side = &(sRxn.products);

			bool bFound = false;
			for (auto &pair: *side)
			{
				if (pair.first == nID)
				{
					pair.second += nCoeff;
					bFound = true;
					break;
				}
			}
			if (!bFound)
				(*side).push_back(std::make_pair(nID, nCoeff));

			// Adding to initial if not defined previously
			bFound = false;
			for (auto &pair : sInput.initial)
			{
				if (pair.first == nID)
				{
					bFound = true;
					break;
				}
			}
			if (!bFound)
				sInput.initial.push_back(std::make_pair(nID, 0));
		}
		else
			std::cout << "Ignoring unrecognised sequence in reaction definition: " << strWord << std::endl;
	}

	// Classify the reaction
	if (nOrder > 3)
	{
		std::cout << "Only reactions up to third order are supported, skipping" << std::endl;
		sRxn.c = -1.0;
	}
	else if (nOrder == 0)
		sRxn.type = RxnType::Influx;
	else if (nOrder == 1)
		sRxn.type = RxnType::Mono;
	else if (nOrder == 2 && sRxn.reactants.size() == 1)
		sRxn.type = RxnType::A2;
	else if (nOrder == 2 && sRxn.reactants.size() == 2)
		sRxn.type = RxnType::AB;
	else if (nOrder == 3 && sRxn.reactants.size() == 1)
		sRxn.type = RxnType::A3;
	else if (nOrder == 3 && sRxn.reactants.size() == 2)
	{
		sRxn.type = RxnType::A2B;
		if (sRxn.reactants[1].second > sRxn.reactants[0].second)
		{
			// Ensure that A is stored first, eliminates need for checking
			auto Temp = sRxn.reactants[0];
			sRxn.reactants[0] = sRxn.reactants[1];
			sRxn.reactants[1] = Temp;
		}
	}
	else if (nOrder == 3 && sRxn.reactants.size() == 3)
		sRxn.type = RxnType::ABC;
	else
		// No idea how we would get here, but just to be sure
		sRxn.c = -1.0;
}

Reaction analyseLine(string strLine, Input &sInput)
{
	Reaction sRxn;

	// Remove comments
	for (unsigned int nChar = 0; nChar < strLine.size(); ++nChar)
	{
		if (strLine[nChar] == '#')
		{
			strLine = strLine.substr(0, nChar);
			break;
		}
	}

	// Split into words
	std::stringstream SStream;
	vector<string> vsWords;
	SStream << strLine;
	while (!SStream.eof())
	{
		string strWord;
		SStream >> strWord;
		// Only non-empty words
		if (!std::all_of(strWord.begin(), strWord.end(), std::isspace))
			vsWords.push_back(strWord);
	}

	if (vsWords.empty())
		// Default Reaction will have negative c and won't be included
		return sRxn;

	// Handle keywords
	if (vsWords[0] == "INIT")
	{
		unsigned int nID = std::atoi(vsWords[1].c_str());
		unsigned int nCount = 1;
		if (vsWords.size() > 2)
			nCount = std::atoi(vsWords[2].c_str());
		bool bFound = false;
		for (auto &pair : sInput.initial)
		{
			if (pair.first == nID)
			{
				pair.second = nCount;
				bFound = true;
				break;
			}
		}
		if (!bFound)
			sInput.initial.push_back(std::make_pair(nID, nCount));
		std::cout << "Starting count of " << vsWords[1] << ": " << vsWords[2] << std::endl;
	}
	else if (vsWords[0] == "DIRECT")
	{
		sInput.first = false;
		std::cout << "Using direct version of the algorithm" << std::endl;
	}
	else if (vsWords[0] == "FIRST")
	{
		sInput.first = true;
		std::cout << "Using first reaction version of the algorithm" << std::endl;
	}
	else if (vsWords[0] == "BINARY")
	{
		sInput.binary = true;
		std::cout << "Using binary search when deciding on a reaction" << std::endl;
	}
	else if (vsWords[0] == "STEPS")
	{
		sInput.steps = std::atoi(vsWords[1].c_str());
		std::cout << "Using " << vsWords[1] << " steps" << std::endl;
	}
	else if (vsWords[0] == "PROGRESS")
	{
		sInput.progress = std::atoi(vsWords[1].c_str());
		std::cout << "Reporting progress every " << vsWords[1] << " steps" << std::endl;
	}
	else if (vsWords[0] == "FPROGRESS")
	{
		sInput.fprogress = std::atoi(vsWords[1].c_str());
		std::cout << "Saving into file every " << vsWords[1] << " steps" << std::endl;
	}
	else if (vsWords[0] == "OUTFILE")
	{
		sInput.outfile = vsWords[1];
		std::cout << "Saving data to file \"" << vsWords[1] << "\"" << std::endl;
	}
	else if (vsWords[0] == "ESTIMATE")
	{
		sInput.estimate = true;
		if (vsWords.size() > 1)
			sInput.maxsteps = std::atoi(vsWords[1].c_str());
		if (vsWords.size() > 2)
			sInput.stepsize = std::atof(vsWords[2].c_str());
		std::cout << "Estimating equilibrium probabilities using at most " << 
			sInput.maxsteps << " steps of size " << sInput.stepsize << std::endl;
	}
	else if (std::isdigit(vsWords[0][0]) || vsWords[0][0] == '-')
		registerReaction(sRxn, vsWords, sInput);
	else
		std::cout << "Ignoring unrecognised sequence: " << strLine << std::endl;

	return sRxn;
}

vector<Reaction> loadInput(string strFilename, Input &sInput)
{
	vector<Reaction> vsReactions;

	std::ifstream InputFile(strFilename.c_str());
	string strLine;
	while (std::getline(InputFile, strLine))
	{
		Reaction sRxn = analyseLine(strLine, sInput);
		if (sRxn.c >= 0.0)
			vsReactions.push_back(sRxn);
	}

	return vsReactions;
}