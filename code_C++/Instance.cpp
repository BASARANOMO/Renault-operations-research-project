#include "Instance.h"
#include <fstream>
#include <iostream>

void instanceError(string const & message)
{
	cerr << "Instance parsing error:" << endl << message << endl;
	throw;
}

void readNameValue(ifstream & read, string letter, int & value, string message){
	string dustbin;
	read >> dustbin >> value;
	if (letter != dustbin) instanceError(message + " erreur: quantité attendue : " + letter \
		+ " quantité obtenue " + dustbin);
}

void readName(ifstream & read, string letter, string message){
	string dustbin;
	read >> dustbin;
	if (letter != dustbin) instanceError(message + " erreur: quantité attendue : " + letter \
		+ " quantité obtenue " + dustbin);
}

void checkValue(int valueRead, int valueAttendue, string message)
{
	if (valueRead != valueAttendue)
	{
		cerr << message << " Erreur:" << endl;
		cerr << "Quantite attendue: " << valueAttendue << "Quantite lue: " << valueRead;
		throw;
	}
}

Instance::Instance(string const & filename, bool verb)
{
	verbose = verb;
	auto read = ifstream(filename);
	if (!read.is_open())
	{
		cerr << "unable to open instance file: " << filename << endl;
	}
	string dustbin = "nothing read up to now";
	readNameValue(read, "J", horizon, "Horizon");
	readNameValue(read, "U", nbUsines, "nbUsines");
	readNameValue(read, "F", nbFournisseurs, "nbFournisseurs");
	readNameValue(read, "E", nbEmballages, "nbEmballages");
	readNameValue(read, "L", metrageLineaireCamion, "metrageLineaireCamion");
	readNameValue(read, "Gamma", coutKilometrique, "coutKilometrique");
	readNameValue(read, "CCam", coutCamion, "coutCamion");
	readNameValue(read, "CStop", coutArretCamion, "coutCamion");

	metrageLineaireEmballages = vector<int>(nbEmballages,-1);
	for (int i = 0; i < nbEmballages; ++i)
	{
		readName(read, "e", "Emballage");
		readName(read, to_string(i), "Emballage");
		readName(read, "l", "Emballage");
		read >> metrageLineaireEmballages[i];
	}

	usines = vector<Usine>(nbUsines);
	for (int i = 0; i < nbUsines; ++i)
	{
		auto & usine = usines[i];
		readNameValue(read, "u", usine.id, "usine id");
		readNameValue(read, "v", usine.vertexId, "usine vertexId");
		readName(read, "coor", "usine latitude et longitude");
		read >> usine.longitude;
		read >> usine.longitude;
		checkValue(usine.id, usine.vertexId, "id et vertexId identiques pour les usines");
		readName(read, "ce", "usine, constantes emballages");
		usine.stocksInitiaux = vector<int>(nbEmballages,-1);
		usine.coutStockExcedentaire = vector<int>(nbEmballages, -1);
		for (int e = 0; e < nbEmballages; ++e)
		{
			readName(read, "e", "Emballage");
			readName(read, to_string(e), "Emballage");
			readNameValue(read, "cr",usine.coutStockExcedentaire[e], "Usine, coutStockExcedentaire");
			readNameValue(read, "b", usine.stocksInitiaux[e], "Usine, stocksInitiaux");
		}
		usine.liberation = vector<vector<int> >(horizon, vector<int>(nbEmballages,-1));
		usine.stockMax = vector<vector<int> >(horizon, vector<int>(nbEmballages,-1));
		readName(read, "lib", "usine, liberation d'emballages ");
		for (int j = 0; j < horizon; ++j)
		{
			readName(read, "j", "usine, demande du jour");
			readName(read, to_string(j), "erreur sur le jour attendu");
			for (int e = 0; e < nbEmballages; ++e)
			{
				readName(read, "e", "Emballage");
				readName(read, to_string(e), "Emballage");
				readNameValue(read, "b", usine.liberation[j][e], "liberation du jour " \
					+ to_string(j) + " pour l'emballage " + to_string(e));
				readNameValue(read, "r", usine.stockMax[j][e], "stockMax du jour " \
					+ to_string(j) + " pour l'emballage " + to_string(e));
			}			
		}
	}

	fournisseurs = vector<Fournisseur>(nbFournisseurs);
	for (int i = 0; i < nbFournisseurs; ++i)
	{
		auto & fournisseur = fournisseurs[i];
		readNameValue(read, "f", fournisseur.id, "fournisseur id");
		readNameValue(read, "v", fournisseur.vertexId, "fournisseur vertexId");
		readName(read, "coor", "latitude et longitude fournisseur");
		read >> fournisseur.latitude;
		read >> fournisseur.longitude;
		checkValue(fournisseur.id + nbUsines, fournisseur.vertexId, "id et vertexId pour les fournisseurs");
		readName(read, "ce", "fournisseur, constantes emballages");
		fournisseur.stocksInitiaux = vector<int>(nbEmballages,-1);
		fournisseur.coutStockExcedentaire = vector<int>(nbEmballages, -1);
		fournisseur.coutExpeditionCarton = vector<int>(nbEmballages, -1);
		for (int e = 0; e < nbEmballages; ++e)
		{
			readName(read, "e", "Emballage");
			readName(read, to_string(e), "Emballage");
			readNameValue(read, "cr",fournisseur.coutStockExcedentaire[e], "fournisseur, coutStockExcedentaire");
			readNameValue(read, "cexc",fournisseur.coutExpeditionCarton[e], "fournisseur, coutStockExcedentaire");
			readNameValue(read, "b", fournisseur.stocksInitiaux[e], "fournisseur, stocksInitiaux");
		}
		fournisseur.demande = vector<vector<int> >(horizon, vector<int>(nbEmballages,-1));
		fournisseur.stockMax = vector<vector<int> >(horizon, vector<int>(nbEmballages,-1));
		readName(read, "dem", "fournisseur, demande d'emballages ");
		for (int j = 0; j < horizon; ++j)
		{
			readName(read, "j", "fournisseur, demande du jour");
			readName(read, to_string(j), "erreur sur le jour attendu");
			for (int e = 0; e < nbEmballages; ++e)
			{
				readName(read, "e", "Emballage");
				readName(read, to_string(e), "Emballage");
				readNameValue(read, "b", fournisseur.demande[j][e], "demande du jour " \
					+ to_string(j) + " pour l'emballage " + to_string(e));
				readNameValue(read, "r", fournisseur.stockMax[j][e], "stockMax du jour " \
					+ to_string(j) + " pour l'emballage " + to_string(e));
			}			
		}
	}
	int nbVertices = nbFournisseurs + nbUsines;
	distanceMatrix = vector<vector<int> >(nbVertices, vector<int>(nbVertices, -1));
	for (int o = 0; o < nbVertices; ++o)
	{
		for (int d = 0; d < nbVertices; ++d)
		{
			readName(read, "a", "arc");
			readName(read, to_string(o), "arc origin");
			readName(read, to_string(d), "arc destination");
			readNameValue(read, "d", distanceMatrix[o][d], "arc distance");
		}
	}	
	read.close();
}