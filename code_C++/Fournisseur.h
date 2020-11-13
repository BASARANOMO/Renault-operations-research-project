#ifndef __FOURNISSEUR__
#define __FOURNISSEUR__

#include <vector>
using namespace std;

class Fournisseur
{
public:
	int id;
	int vertexId;
	double latitude;
	double longitude;
	vector<int> stocksInitiaux;
	vector<int> coutStockExcedentaire;
	vector<int> coutExpeditionCarton;
	vector<vector<int> > demande;			// demande[jour][emballage]
	vector<vector<int> > stockMax;	

	Fournisseur();
	~Fournisseur() = default;

	void print();
	
};

#endif