#include <vector>
#include <algorithm>
#include <iostream>
#include "Integrate.hpp"
#include <random>


using namespace std;

Integrate::Integrate(double (*func)(double, double, double), std::vector< double > ranges, int MCSteps){
	fun_ = func;
	_ranges = ranges;
	_MCSteps = MCSteps;
	_dim = _ranges.size()/2;
	/*
	 * La questione dei Bin in X e Z va scelta in modo da mantenere
	 * la densita` di punti montecarlo, per volumetto, costante e ad un livello
	 * necessario per garantire la precisione con cui si vuole stimare l'integrale
	 *
	 */



}


Integrate::~Integrate(){

};

double Integrate::cubatura_naif(){
	
	// x,y sono estratti all'interno di ogni bin
	double x,y,z,w;

	// il volume totale viene costruito aggiungendo di volta in volta i volumetti dei singoli bin.
	double volume = 0;

	// Numero di hit (ho insierito una distizione + e -, ma forse e` ridontante... cercare di
	// capire qual e` l'implementazione migliore)
	int hit=0;
	
	
	// definisco i due generatori delle x e y (MersenneTwister)
	std::random_device rd_x;  
	std::mt19937 gen_x(rd_x()); 
	
	std::random_device rd_y;  
	std::mt19937 gen_y(rd_y());

	std::random_device rd_z;
	std::mt19937 gen_z(rd_z());

	std::random_device rd_w;
	std::mt19937 gen_w(rd_w());
	
	// Il numero di suddivisioni del range X e del range Y di integrazione. 
	// Questo numeor ora e` settato a caso, va ottimizzato in modo da mantenere
	// la densita` di di punti montecarlo per volumetto costante.
	int _N_x_bins = 200;
	int _N_y_bins = 200;
	int _N_z_bins = 200;

	//Passo con cui mi sposto tra i volumetti
	double dx = (_ranges[1] - _ranges[0]) / _N_x_bins;
	double dy = (_ranges[3] - _ranges[2]) / _N_y_bins;
	double dz = (_ranges[5] - _ranges[4]) / _N_z_bins;


	for(int i=0; i < _N_x_bins; i++){
		double x_i = _ranges[0]+i*dx;
		double x_f = _ranges[0]+(i+1)*dx;
		for(int l=0; l< _N_z_bins; l++){
			double z_i = _ranges[4]+l*dz;
			double z_f = _ranges[4]+(l+1)*dz;
			for(int k=0; k < _N_y_bins; k++){
				double y_i = _ranges[2]+k*dy;
				double y_f = _ranges[2]+(k+1)*dy;
			
				// Una volta che ho x_i, x_f, y_i, y_f posso definire i ranges
				// delle distribuzioni di x e y			
				std::uniform_real_distribution<> dis_x(x_i, x_f); 
				std::uniform_real_distribution<> dis_y(y_i, y_f);
				std::uniform_real_distribution<> dis_z(z_i, z_f);
			
				//Va definito meglio il range di z, ora fa schifo
				std::vector <double> w_v;
				w_v.push_back(fun_(x_i,y_i,z_i));
				w_v.push_back(fun_(x_f,y_f,z_f));
				w_v.push_back(fun_(x_i+dx/2.,y_i+dy/2.,z_f+dz/2.));
				w_v.push_back(fun_(x_i,y_f,z_i));
				w_v.push_back(fun_(x_f,y_i,z_f));

				//z_v.push_back(fun_((x_f-x_i)/2.,(y_f-y_i)/2.));
				double w_i = w_v[std::distance(w_v.begin(),std::max_element(std::begin(w_v), std::end(w_v)))];
				double w_f = w_v[std::distance(w_v.begin(),std::min_element(std::begin(w_v), std::end(w_v)))];
				double inf = w_i/10.;
				double sup = w_f/10.;
				std::uniform_real_distribution<> dis_w(w_i-inf, w_f+sup); 
				w_v.erase(w_v.begin(), w_v.end());
	
				// per ogni volumetto, il numero degli hit si deve azzerare.
				hit=0;
				for(int smc=0; smc < _MCSteps; smc++){
					x = dis_x(gen_x);
					y = dis_y(gen_y);
					z = dis_z(gen_z);
					w = dis_w(gen_w);
					if( w < fun_(x,y,z) ) hit++;
				}
				volume += dx*dy*dz*((w_f+sup)-(w_i-inf))*(double(double(hit)/double(_MCSteps))); 
				volume += dx*dy*dz*(w_i-inf);
			}
		}
		//cout << "Completed percentage: " << double(i)/double(_N_x_bins) << endl;
	}

	return volume;

}




