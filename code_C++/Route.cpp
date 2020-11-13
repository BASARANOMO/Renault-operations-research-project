#include "Route.h"
#include "cassert"
#include <sstream>
#include <iostream>


Route::Route()
{
	id = "not_initialized";
	usine = -1;
}

// void routeError(string const & r){
// 	cerr << "Error reading route:" << endl;
// 	cerr << r << endl;
// 	throw;
// }

// Route::Route(string const & routeLine)
// {
// 	auto read = stringstream(routeLine);
// 	string dustbin;
// 	read >> dustbin;
// 	if (dustbin != "r") routeError(routeLine);
// 	read >> id;
// 	read >> dustbin;
// 	if (dustbin != "j") routeError(routeLine);
// 	read >> jour >> dustbin;
// 	if (dustbin != "x") routeError(routeLine);
// 	read >> nbRealisations >> dustbin;
// 	if (dustbin != "F") routeError(routeLine);
// 	read >> nbFournisseurs;
// 	if (nbFournisseurs < 1) routeError(routeLine + " route avec trop peu de founisseurs");
// 	if (nbFournisseurs > 4) routeError(routeLine + " route avec trop de fournisseurs");
	
	

// }