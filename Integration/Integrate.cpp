#include <vector>
#include <algorithm>
#include <iostream>
#include "Integrate.hpp"
#include <random>


using namespace std;

Integrate::Integrate(double (*func)(double, double), std::vector< double > ranges, int MCSteps){
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
	double x,y,z;

	// il volume totale viene costruito aggiungendo di volta in volta i volumetti dei singoli bin.
	double volume = 0;

	// Numero di hit (ho insierito una distizione + e -, ma forse e` ridontante... cercare di
	// capire qual e` l'implementazione migliore)
	int hit_plus=0, hit_minus=0;
	
	
	// definisco i due generatori delle x e y (MersenneTwister)
	std::random_device rd_x;  
	std::mt19937 gen_x(rd_x()); 
	
	std::random_device rd_y;  
	std::mt19937 gen_y(rd_y());

	std::random_device rd_z;
	std::mt19937 gen_z(rd_z());
	
	// Il numero di suddivisioni del range X e del range Y di integrazione. 
	// Questo numeor ora e` settato a caso, va ottimizzato in modo da mantenere
	// la densita` di di punti montecarlo per volumetto costante.
	int _N_x_bins = 100;
	int _N_y_bins = 100;

	//Passo con cui mi sposto tra i volumetti
	double dx = (_ranges[1] - _ranges[0]) / _N_x_bins;
	double dy = (_ranges[3] - _ranges[2]) / _N_y_bins;


	for(int i=0; i < _N_x_bins; i++){
		double x_i = _ranges[0]+i*dx;
		double x_f = _ranges[0]+(i+1)*dx;
		for(int k=0; k < _N_y_bins; k++){
			double y_i = _ranges[2]+k*dy;
			double y_f = _ranges[2]+(k+1)*dy;
			
			// Una volta che ho x_i, x_f, y_i, y_f posso definire i ranges
			// delle distribuzioni di x e y			
			std::uniform_real_distribution<> dis_x(x_i, x_f); 
			std::uniform_real_distribution<> dis_y(y_i, y_f);
			
			//Va definito meglio il range di z, ora fa schifo
			std::vector <double> z_v;
			z_v.push_back(fun_(x_i,y_i));
			z_v.push_back(fun_(x_f,y_f));
			double z_i = z_v[std::distance(z_v.begin(),std::max_element(std::begin(z_v), std::end(z_v)))];
			double z_f = z_v[std::distance(z_v.begin(),std::min_element(std::begin(z_v), std::end(z_v)))];
			std::uniform_real_distribution<> dis_z(z_i-(z_i/10.), z_f+(z_f/10.)); //Bho, lo allungo del 10%...
			z_v.erase(z_v.begin(), z_v.end());
	
			// per ogni volumetto, il numero degli hit si deve azzerare.
			hit_plus=0;
			for(int smc=0; smc < _MCSteps; smc++){
				x = dis_x(gen_x);
				y = dis_y(gen_y);
				z = dis_z(gen_z);
				if( z < fun_(x,y) ) hit_plus++;
			}
			volume += (x_f-x_i)*(y_f-y_i)*((z_f+(z_f/10.))-(z_i-(z_i/10.)))*(double(double(hit_plus)/double(_MCSteps))); 
			volume += (x_f-x_i)*(y_f-y_i)*(z_i-(z_i/10.));
		}
	}

	return volume;

}




