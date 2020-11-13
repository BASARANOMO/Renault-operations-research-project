#ifndef __USINE__
#define __USINE__
#include <vector>

using namespace std;

class Usine
{
public:
	int id;
	int vertexId;
	double latitude;
	double longitude;
	vector<int> stocksInitiaux;
	vector<int> coutStockExcedentaire;
	vector<vector<int> > liberation;	// liberation[jour][emballage]
	vector<vector<int> > stockMax;

	Usine();
	~Usine() = default;
	
};

#endif
