#include <vector>
#include <iostream>
#include "Integrate.hpp"
#include <random>


using namespace std;

Integrate::Integrate(double (*func)(double, double), std::vector< double > ranges, int MCSteps){
	fun_ = func;
	_ranges = ranges;
	_MCSteps = MCSteps;

	//usare qualche metodo di ranges per ottenere la dim
	
	_dim = _ranges.size()/2;

	cout << _dim << endl;

}


Integrate::~Integrate(){

};




double Integrate::cubatura_naif(){

	double x,y,z;
	double z_f;
	int hit_plus=0, hit_minus=0;

	std::random_device rd_x;  
	std::mt19937 gen_x(rd_x()); 
	std::uniform_real_distribution<> dis_x(_ranges[0], _ranges[1]);
	
	std::random_device rd_y;  
	std::mt19937 gen_y(rd_y()); 
	std::uniform_real_distribution<> dis_y(_ranges[2], _ranges[3]);
	
	std::random_device rd_z;  
	std::mt19937 gen_z(rd_z()); 
	std::uniform_real_distribution<> dis_z(_ranges[4], _ranges[5]);

	for(int i=0; i<_MCSteps; i++){
		x = dis_x(gen_x);
		y = dis_y(gen_y);
		z = dis_z(gen_z);
		z_f = fun_(x,y);

		if( z <= z_f ) hit_plus++;
	}

	return ((_ranges[1]-_ranges[0])*(_ranges[3]-_ranges[2])*(_ranges[5]-_ranges[4]) ) * double(double(hit_plus)/double(_MCSteps))  ;

}




