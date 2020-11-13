#ifndef __INSTANCE__
#define __INSTANCE__

#include "Usine.h"
#include "Fournisseur.h"

#include <string>

class Instance
{
public:
	bool verbose;

	int horizon;
	int nbUsines;
	int nbFournisseurs;
	int nbEmballages;
	int metrageLineaireCamion;
	int coutKilometrique;
	int coutCamion;
	int coutArretCamion;
	
	vector<int> metrageLineaireEmballages;
	vector<Usine> usines;
	vector<Fournisseur> fournisseurs;
	vector<vector<int> > distanceMatrix;

	Instance(string const &, bool);
	~Instance() = default;
	
};

#endif