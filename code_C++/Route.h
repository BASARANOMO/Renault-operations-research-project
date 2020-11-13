#ifndef __ROUTE__
#define __ROUTE__

#include <string>
#include <vector>

using namespace std;



class Route
{
public:
	string id;
	int jour;
	int nbRealisations;
	int usine;
	int nbFournisseurs;
	vector<int> fournisseurs;
	vector<vector<int> > quantitesLivrees; 	// quantitesLivrees[fournisseur][emballage] 

	Route();
	Route(string const &);
	~Route() = default;
};

#endif