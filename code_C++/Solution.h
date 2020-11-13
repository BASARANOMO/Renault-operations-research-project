#ifndef __SOLUTION__
#define __SOLUTION__

#include "Instance.h"
#include "Route.h"

class Solution
{
public:
	Instance const & instance;

	// Read in file
	int nbRoutes;
	vector<Route> routes;

	// Stocks recomputed
	vector<vector<vector<int> > > stockUsineLeSoir;			// [usine][jour][emballage]
	vector<vector<vector<int> > > stockFournisseurLeSoir;	// [fournisseur][jour][emballage]
	int cost;

	Solution(Instance const & inst, string const & solutionFileName);
	~Solution() = default;
	
	bool calculDesCoutsEtVerificationsFaisabilite();

private:
	bool calculStocksEtCoutsUsines();
	bool calculStocksEtCoutsFournisseurs();
	bool calculCoutsRoutes();

};

#endif