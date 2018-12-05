#include<vector>
#include <stdlib.h>
#include<iostream>
#include<random>


#include"xorshift.h"

using namespace std;

const int time_obs = 1000000;
const int N_det    = 10;


void printTOD(std::vector< std::vector< double > > TOD){
	for(auto i : TOD){
		for(auto k : i)
			cout<< ' ' << k;
		cout<<endl;
	}
}


int main(){

	// Define TOD
	std::vector< std::vector< double > > TOD;
	
	// Define gaussian pseudo-random number generator 
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(5.0,2.0);

	// Start observation
	for(int t=0; t < time_obs; t++){
		std::vector<double> v_p;
		for(int p =0; p < N_det; p++){
			v_p.push_back(distribution(generator));
		}
		TOD.push_back(v_p);
	}
	printTOD(TOD);
	return 0;
}







