#include "Solution.h"
#include <iostream>
#include <fstream>
#include <algorithm>

void solutionError(string const & message)
{
	cerr << "Solution parsing error:" << endl << message << endl;
	throw;
}

void readSolutionNameValue(ifstream & read, string letter, int & value, string message){
	string dustbin;
	read >> dustbin >> value;
	if (letter != dustbin) solutionError(message + " erreur: quantité attendue : " + letter \
		+ " quantité obtenue " + dustbin);
}

void readSolutionName(ifstream & read, string letter, string message){
	string dustbin;
	read >> dustbin;
	if (letter != dustbin) solutionError(message + " erreur: quantité attendue : " + letter \
		+ " quantité obtenue " + dustbin);
}

Solution::Solution(Instance const & inst, string const & solutionFileName) : instance(inst)
{
	auto read = ifstream(solutionFileName);
	if (!read.is_open()) solutionError("Unable to open " + solutionFileName);
	readSolutionNameValue(read, "R", nbRoutes, "nbRoutes");
	routes = vector<Route>(nbRoutes);
	for (int r = 0; r < nbRoutes; ++r)
	{
		auto & route = routes[r];
		readSolutionName(read, "r", "Route");
		read >> route.id;
		readSolutionNameValue(read, "j", route.jour, "jour de la route " + route.id);
		readSolutionNameValue(read, "x", route.nbRealisations, "variable x indiquant le nombre de fois que la route "\
			+ route.id + " est prise");
		readSolutionNameValue(read, "u", route.usine, "usine de la route " + route.id);
		readSolutionNameValue(read, "F", route.nbFournisseurs, "nombres de fournisseurs sur la route " + route.id);
		route.fournisseurs = vector<int>(route.nbFournisseurs, -1);
		route.quantitesLivrees = vector<vector<int> >(route.nbFournisseurs, vector<int>(instance.nbEmballages, -1));
		for (int f = 0; f < route.nbFournisseurs; ++f)
		{
			readSolutionNameValue(read, "f", route.fournisseurs[f], "Fournisseur numéro " + to_string(f) + " de la route " + route.id);
			for (int e = 0; e < instance.nbEmballages; ++e)
			{
				readSolutionName(read, "e", "Emballage " + to_string(e) + " du fournisseur numéro " + to_string(f)\
					+ "de la route " + route.id);
				readSolutionName(read, to_string(e), "Emballage " + to_string(e) + " du fournisseur numéro " + to_string(f)\
					+ "de la route " + route.id);
				readSolutionNameValue(read, "q", route.quantitesLivrees[f][e], "quantité d'emballage " + to_string(e) +\
					+ "du fournisseur " + to_string(f) + " de la route " + route.id);
			}
		}
	}
}

bool Solution::calculDesCoutsEtVerificationsFaisabilite()
{
	sort(routes.begin(), routes.end(),[](Route const & r1, Route const & r2)
		{
			return r1.jour < r2.jour;
		});
	cost = 0;
	if (!calculStocksEtCoutsUsines()) {
		cerr << "Solution Usine non réalisable" << endl; 	
		return false;
	}
	if (!calculStocksEtCoutsFournisseurs())
	{
		cerr << "Solution fournisseur non réalisable" << endl;
		return false;
	}
	if (!calculCoutsRoutes())
	{
		cerr << "Route non réalisable" << endl;
		return false;
	}
	return true;
}

bool Solution::calculStocksEtCoutsUsines()
{
	stockUsineLeSoir = vector<vector<vector<int> > > (instance.nbUsines,\
		vector<vector<int> >(instance.horizon,\
			vector<int>(instance.nbEmballages,-1)));
	int j = 0;
	int r = 0;
	int costUsine = 0.0;
	// Stocks le matin du premier jour
	for (int u = 0; u < instance.nbUsines; ++u)
	{
		for (int e = 0; e < instance.nbEmballages; ++e)
		{
			stockUsineLeSoir[u][j][e] = instance.usines[u].stocksInitiaux[e];
		}
	}
	while (j < instance.horizon)
	{
		// Libération d'emballages
		for (int u = 0; u < instance.nbUsines; ++u)
		{
			for (int e = 0; e < instance.nbEmballages; ++e)
			{
				stockUsineLeSoir[u][j][e] += instance.usines[u].liberation[j][e];
			}
		}		
		// Envoi sur les routes
		while (r < routes.size() && routes[r].jour == j)
		{
			auto const & route = routes[r];
			for (int f = 0; f < route.nbFournisseurs; ++f)
			{
				for (int e =0; e < instance.nbEmballages; ++e)
				{
					stockUsineLeSoir[route.usine][j][e] -= \
						route.nbRealisations *route.quantitesLivrees[f][e];
				}
			}
			++r;
		}
		// Faisabilite et couts des stocks
		for (int u = 0; u < instance.nbUsines; ++u)
		{
			for (int e = 0; e < instance.nbEmballages; ++e)
			{
				if (stockUsineLeSoir[u][j][e] < 0)
				{
					cerr << "L'usine " + to_string(u) + " a un stock négatif d'emballage " \
						+ to_string(e) + " le soir du jour " + to_string(j) << endl;
					return false;
				}
				if (stockUsineLeSoir[u][j][e] > instance.usines[u].stockMax[j][e])
				{
					costUsine += instance.usines[u].coutStockExcedentaire[e] * \
						(stockUsineLeSoir[u][j][e] - instance.usines[u].stockMax[j][e]);
				}
			}
		}		
		++j;
		// Stocks le matin du jour suivant
		if (j < instance.horizon)
		{
			for (int u = 0; u < instance.nbUsines; ++u)
			{
				for (int e = 0; e < instance.nbEmballages; ++e)
				{
					stockUsineLeSoir[u][j][e] = stockUsineLeSoir[u][j-1][e];
				}
			}			
		}
	}
	if (r != routes.size())
	{
		cerr << "Toutes les routes ne partent pas un jour de l'horizon" << endl;
		return false;
	}
	if (instance.verbose) cout << "cout usines: " << costUsine << endl;
	cost += costUsine;
	return true;
}

bool Solution::calculStocksEtCoutsFournisseurs()
{
	stockFournisseurLeSoir = vector<vector<vector<int> > > (instance.nbFournisseurs,\
		vector<vector<int> >(instance.horizon,\
			vector<int>(instance.nbEmballages,-1)));
	int j = 0;
	int r = 0;
	int costFournisseur = 0;
	vector<vector<int> > costFournisseurJour = vector<vector<int> >(instance.nbFournisseurs, vector<int>(instance.horizon, 0));

	for (int f = 0; f < instance.nbFournisseurs; ++f)
	{
		for (int e = 0; e < instance.nbEmballages; ++e)
		{
			stockFournisseurLeSoir[f][j][e] = instance.fournisseurs[f].stocksInitiaux[e];
		}
	}
	while (j < instance.horizon)
	{
		// Consommation d'emballages et cost expedition carton
		for (int f = 0; f < instance.nbFournisseurs; ++f)
		{
			for (int e = 0; e < instance.nbEmballages; ++e)
			{
				stockFournisseurLeSoir[f][j][e] -= instance.fournisseurs[f].demande[j][e];
				if (stockFournisseurLeSoir[f][j][e] < 0)
				{
					costFournisseur += instance.fournisseurs[f].coutExpeditionCarton[e] * \
						(- stockFournisseurLeSoir[f][j][e]);
					costFournisseurJour[f][j] += instance.fournisseurs[f].coutExpeditionCarton[e] * \
						(- stockFournisseurLeSoir[f][j][e]);
					stockFournisseurLeSoir[f][j][e] = 0;
				}
 			}
		}

		// Reception depuis les routes
		while (r < routes.size() && routes[r].jour == j)
		{
			auto const & route = routes[r];
			for (int f_r = 0; f_r < route.nbFournisseurs; ++f_r)
			{
				for (int e =0; e < instance.nbEmballages; ++e)
				{
					if (route.quantitesLivrees[f_r][e] < 0)
					{
						cerr << "La route " << route.id << "livre une quantité négative \
							d'emballage à son fournisseur numéro " << to_string(f_r) << \
							" le jour " << to_string(j) << endl; 
						return false;
					}
					stockFournisseurLeSoir[route.fournisseurs[f_r]][j][e] += \
						route.nbRealisations * route.quantitesLivrees[f_r][e];
				}
			}
			++r;
		}
		// Faisabilite et couts des stocks
		for (int f = 0; f < instance.nbFournisseurs; ++f)
		{
			for (int e = 0; e < instance.nbEmballages; ++e)
			{
				if (stockFournisseurLeSoir[f][j][e] < 0)
				{
					cerr << "Le fournisseur " + to_string(f) + " a un stock négatif d'emballage " \
						+ to_string(e) + " le soir du jour " + to_string(j) << endl;
					return false;
				}
				if (stockFournisseurLeSoir[f][j][e] > instance.fournisseurs[f].stockMax[j][e])
				{
					costFournisseur += instance.fournisseurs[f].coutStockExcedentaire[e] * \
						(stockFournisseurLeSoir[f][j][e] - instance.fournisseurs[f].stockMax[j][e]);
					costFournisseurJour[f][j] += instance.fournisseurs[f].coutStockExcedentaire[e] * \
						(stockFournisseurLeSoir[f][j][e] - instance.fournisseurs[f].stockMax[j][e]);
				}
			}
		}		
		++j;
		// Stock le matin du jour suivant
		if (j < instance.horizon)
		{
			for (int f = 0; f < instance.nbFournisseurs; ++f)
			{
				for (int e = 0; e < instance.nbEmballages; ++e)
				{
					stockFournisseurLeSoir[f][j][e] = stockFournisseurLeSoir[f][j-1][e];
				}
			}		
		}
	}
	if (r != routes.size())
	{
		cerr << "Toutes les routes ne partent pas un jour de l'horizon" << endl;
		return false;
	}
	if (instance.verbose) cout << "cout fournisseurs : " << costFournisseur << " avec le détail :" << endl;
	for (int f = 0; f < instance.nbFournisseurs; f++)
	{
		for (int j =0; j< instance.horizon; ++j){
			if (instance.verbose) cout << "  Fournisseur " << f << " jour " << j << " cout " << costFournisseurJour[f][j] <<  endl;
		}
	}

	if (instance.verbose) cout << "Détail des stocks le soir: " << endl;
	for (int f = 0; f < instance.nbFournisseurs; f++)
	{
		for (int j =0; j< instance.horizon; ++j){
			for (int e = 0; e< instance.nbEmballages; ++e)
			{
				if (instance.verbose) cout << "  fournisseur " << f << " jour " << j << " emballage " << e << " st " << stockFournisseurLeSoir[f][j][e] << endl; 
			}
		}
	}
	
	cost += costFournisseur;
	return true;
}

bool Solution::calculCoutsRoutes()
{
	int costRoute = 0;
	for (auto && route : routes)
	{
		if (route.nbFournisseurs == 0)
		{
			cerr << "La route " << route.id << " ne contient pas de fournisseurs" << endl;
			return false;
		}
		if (route.nbRealisations == 0)
		{
			cerr << "Warning: La route " << route.id << " n'est par réalisée" << endl;
			continue;
		}
		int vertexUsine = instance.usines[route.usine].vertexId;
		int vertexPremierFourniseurs = instance.fournisseurs[route.fournisseurs[0]].vertexId;
		int distance = instance.distanceMatrix[vertexUsine][vertexPremierFourniseurs];

		for (int f = 0; f < route.nbFournisseurs - 1; ++f)
		{
			int origine = instance.fournisseurs[route.fournisseurs[f]].vertexId;
			int destination = instance.fournisseurs[route.fournisseurs[f + 1]].vertexId;
			distance += instance.distanceMatrix[origine][destination];
		}

		costRoute += (instance.coutCamion + instance.coutArretCamion * route.nbFournisseurs \
			+ instance.coutKilometrique * distance) * route.nbRealisations;
	}
	if (instance.verbose) cout << "couts routes " << costRoute << endl;
	cost += costRoute;
	return true;
}
