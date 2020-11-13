#include "Instance.h"
#include "Solution.h"
#include <iostream>
#include <limits>
#include <fstream>

string executable = "";

int help(string message)
{
	cerr << message << endl;
	cout << "Usage: " << executable << " instanceName <options>" << endl;
	cout << endl;
	cout << "File instanceName will be read and its format checked. It should \
	 contain an instance of the dispatch problem" << endl;
	cout << "Options: " << endl;
	cout << " -checkSol <string>: read and check the solution provided" << endl;
	cout << " -out <string>: output file name (default = out.txt)" << endl;
	cout << " -v: verbose " << endl;
	return 1;
}

int main(int argc, char const *argv[])
{
	executable = string(argv[0]);
	bool verbose = false;
	if (argc < 2) return help("Too few arguments");
	
	// Parse arguments
	auto instanceFile = string(argv[1]);
	bool checkSolution = false;
	string solutionFile = "notSet";
	string outputFile = "out.txt";
	int i = 2;
	while (i < argc)
	{
		auto option = string(argv[i]);
		if (option == "-checkSol"){
			solutionFile = string(argv[++i]);
			checkSolution = true;
		}
		else if (option == "-out"){
			outputFile = string(argv[++i]);
		}
		else if (option == "-v") verbose = true;
		else return help("Bad option " + option);
		++i;
	}

	auto instance = Instance(instanceFile, verbose);
	cout << "Instance lue. Le format de fichier est respecté" << endl;

	if (checkSolution)
	{
		int cost = numeric_limits<int>::max();
		auto solution = Solution(instance, solutionFile);
		cout << "Solution lue, le format de fichier est respecté" << endl;
		if (!solution.calculDesCoutsEtVerificationsFaisabilite())
		{
			cerr << "La solution n'est pas réalisable" << endl;
		}
		else
		{
			cost = solution.cost;
			cout << "La solution est réalisable" << endl;
		}
		cout << "Le cout de la solution est " << cost << endl;
	
		auto outfile = ofstream(outputFile, std::ios_base::app);
		outfile << instanceFile << " " << solutionFile << " " << cost << endl;
		outfile.close();
	}

	return 0;
}