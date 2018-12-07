#include <stdlib.h>
#include <vector>
#include <iostream>
#include <math.h>
#include "Integrate.hpp"
#include "funct.hpp"

using namespace std;





int main(int argv, char **argc){

	double integral = 0;

	//std::vector<double> ranges = {0,2.*M_PI/3.,0,2.*M_PI/3.,0,2.*M_PI/3.};
	int MCSteps = atoi(argc[1]);
	std::vector<double> ranges = {-10.,10.,-10.,10.,-10.,10.};

	
	//Integrate integ(square, ranges, MCSteps);
	//Integrate integ(sinxy, ranges, MCSteps);
	//Integrate integ(plane_o, ranges, MCSteps);
	Integrate integ(gaussian, ranges, MCSteps);


	integral = integ.cubatura_naif();

	cout << integral << endl;


	return 0;

}




