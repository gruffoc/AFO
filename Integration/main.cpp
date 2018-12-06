#include <stdlib.h>
#include <vector>
#include <iostream>
#include "Integrate.hpp"
#include "funct.hpp"

using namespace std;





int main(){

	std::vector<double> ranges = {0,3,0,3,0,3};
	int MCSteps = 100000000;
	
	

	Integrate integ(square, ranges, MCSteps);

	cout << integ.cubatura_naif() << endl;;



	return 0;

}




