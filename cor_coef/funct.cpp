#include <stdlib.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <functional>

#include <iostream>
#include <algorithm>
#include <iterator>



double gaussian (double *k, size_t dim, void *params){

  (void)(dim); 
  (void)(params);
  return exp(-k[0]*k[0] - k[1]*k[1] -k[2]*k[2]);

}


double omega(double *k){
	// theta_b = 0.17deg - STRIP telescope W-band
	
	double lambda  = 0.002;
	double omega_0 = lambda / (M_PI * 0.002967059 );
	double w;

	std::vector <double> r_s = { cos(k[1])*cos(k[2]), cos(k[1])*sin(k[2]), sin(k[1])    };
	std::vector <double> r   = { k[0]*cos(k[1])*cos(k[2]), k[0]*cos(k[1])*sin(k[2]), k[0]*sin(k[1])   };

	w = omega_0 * omega_0 * (1  +    pow(( (lambda*     std::inner_product(r_s.begin(), r_s.end(), r.begin(),0)    )/(M_PI*omega_0*omega_0)    ),2)    );


	return w; 
}



double Beam (double *k){
	// Convenzione sistema di coordinate
	// k[0] = r; k[1] = theta; k[2] = phi
	
	double lambda = 0.002;
	std::vector <double> r_s = { cos(k[1])*cos(k[2]), cos(k[1])*sin(k[2]), sin(k[1])    };
	std::vector <double> r   = { k[0]*cos(k[1])*cos(k[2]), k[0]*cos(k[1])*sin(k[2]), k[0]*sin(k[1])   };

	double part_1 = (2 * pow(lambda,2.) * pow(std::inner_product(r_s.begin(), r_s.end(), r.begin(), 0),2) ) / (M_PI * pow(omega(k),2));
	double part_2 = expf(-( 2*(std::inner_product(r.begin(), r.end(), r.begin(),0) - std::inner_product(r_s.begin(), r_s.end(), r.begin(), 0)  ) / (pow(omega(k),2))));
	
	//std::cout << "Ciao" << std::endl;
	//std::cout << part_1 << ' ' << part_2 << ' ' << part_1*part_2 << std::endl;
	return part_1*part_2;
}


double Chi_1(double *k, double *k2){

	double X1_0=1.;
	double L_0 = 300;
	std::vector <double> r1   = { k[0]*cos(k[1])*cos(k[2]), k[0]*cos(k[1])*sin(k[2]), k[0]*sin(k[1])   };
	std::vector <double> r2   = { k2[0]*cos(k2[1])*cos(k2[2]), k2[0]*cos(k2[1])*sin(k2[2]), k2[0]*sin(k2[1])   };

	std::vector <double> diff;

	std::set_difference(r1.begin(), r1.end(), r2.begin(), r2.end(), std::back_inserter(diff) );


	
	return X1_0 * expf( -  std::inner_product(diff.begin(), diff.end(), diff.begin(),0)  / (2*L_0*L_0) );

}

double Chi_2(double *k, double *k2){
	double X2_0 = 1.;
	double z_0  = 2500;
	double z, z_prime;

	z = k[0]*sin(k[1]);
	z_prime = k2[0]*sin(k2[1]);

	return X2_0 * expf( - ( z+z_prime )/(2*z_0)   );

}

double T_p(double *k){
	
	double T_g = 280;
	double z_atm = 40000;
	double z;

	z = k[0]*sin(k[1]);

	return T_g*(1-(z/z_atm));

}

double Jac(double *k){
	
	std::vector <double> r   = { k[0]*cos(k[1])*cos(k[2]), k[0]*cos(k[1])*sin(k[2]), k[0]*sin(k[1])   };

	return pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],0.5)*cos(k[1]);
}



double C00_t (double *k, size_t dim, void *time){
	double t=0;
	double lambda = 0.002;
	double f = 0.01; //in Hz
	double v1[3], v2[3];

	v1[0] = k[0];
	v1[1] = k[1];
	v1[2] = k[2];

	v2[0] = k[3];
	v2[1] = k[4];
	v2[2] = k[5];

	return ( 1. / pow(lambda,4) ) * ( Beam(v1) * Beam(v2) * Chi_1(v1,v2) * Chi_2(v1,v2) * T_p(v1) * T_p(v2) * Jac(v1) * Jac(v2) ); 
}

