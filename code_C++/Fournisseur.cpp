#include "Fournisseur.h"
#include <iostream>

using namespace std;

Fournisseur::Fournisseur()
{
	id = -1;
	vertexId = -1;
}

void Fournisseur::print()
{
	cout << "Fournisseur " << id << endl;
	cout << "  cout stock ";
	for (auto && c : coutStockExcedentaire) cout << " " << c;
	cout << endl;
	cout << " stock initial ";
	for (auto && s : stocksInitiaux) cout << " " << s;
	cout << endl;
	cout << " stock max ";
	for (auto && jour : stockMax)
	{
		for (auto && emb : jour) cout << " " << emb;
		cout << " ";
	}
	cout << endl;

}