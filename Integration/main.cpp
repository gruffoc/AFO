#include <stdlib.h>
#include <vector>
#include <iostream>
#include "Integrate.hpp"
#include "funct.hpp"

using namespace std;





int main(int argv, char **argc){

	std::vector<double> ranges = {0,2,0,2,0,3};
	int MCSteps = atoi(argc[1]);
	
	

	Integrate integ(square, ranges, MCSteps);
	//Integrate integ(sinxy, ranges, MCSteps);
	//Integrate integ(plane_o, ranges, MCSteps);
	//Integrate integ(gaussian, ranges, MCSteps);
	cout << integ.cubatura_naif() << endl;;



	return 0;

}




