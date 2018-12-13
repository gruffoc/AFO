#include <stdlib.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <functional>

double gaussian (double *k, size_t dim, void *params){

  (void)(dim); 
  (void)(params);
  return exp(-k[0]*k[0] - k[1]*k[1] -k[2]*k[2]);

}


double omega(doube *k, size_t dim, void *params){
	double omega_0 = 1.;
	double lambda  = 1.;

	std::vector <double> r_s = { cos(k[1])*cos(k[2]), cos(k[1])*sin(k[2]), sin(k[1])    };
	std::vector <double> r   = { k[0]*cos(k[1])*cos(k[2]), k[0]*cos(k[1])*sin(k[2]), k[0]*sin(k[1])   };

	w = omega_0 * omega_0 * (1  +    pow(( (lambda*     std::inner_product(r_s.begin(), r_s.end(), r.begin(),0)    )/(M_PI*omega_0*omega_0)    ),2)    );


	return w; 
}



double Beam (double *k, size_t dim, void *params){
	// Convenzione sistema di coordinate
	// k[0] = r; k[1] = theta; k[2] = phi
	
	double lambda = 1.;
	std::vector <double> r_s = { cos(k[1])*cos(k[2]), cos(k[1])*sin(k[2]), sin(k[1])    };
	std::vector <double> r   = { k[0]*cos(k[1])*cos(k[2]), k[0]*cos(k[1])*sin(k[2]), k[0]*sin(k[1])   };

	double part_1 = (2 * pow(lamda,2.) * pow(std::inner_product(r_s.begin(), r_s.end(), r.begin(), 0),2) ) / (M_PI * pow(omega(k,3,0),2));
	double part_2 = expf( 2*(std::inner_product(r.begin(), r.end(), r.begin(),0) - std::inner_product(r_s.begin(), r_s.end(), r.begin(), 0)  ) / (pow(omega(k,3,0),2)));

	return part_1*part_2;
}

double Chi_1 (double *k, size_t dim, void *params){

	double chi_10 = 1.;
	double L0     = 1.;

	// questa struttura e` per ora pensata per prendere solo 1 punto spaziale 
	// bisogna ripensarla un po' in modo da inserire 2 punti saziali a cui
	// corrispondono 2 punti temporali

	return 1.0;

}

double Chi_2 (double *k, size_t dim, void *params){

	return 1.0;

}

double T_p (double *k, size_t dim, void *params){

	return 1.0;

}

double jacobian(double *k, size_t dim, void *params){

	
}


double C00_t (double *k, size_t dim, void *params){


	return 1.0;
}

