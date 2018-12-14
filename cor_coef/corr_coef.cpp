#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <iostream>
#include "funct.hpp"




int main (int argc, char* argv[]){

  double res, err;


  // Integration ranges xl -> lower and xu -> upper
  // La convenzione di coordinate dovrebbe essere: R \theta e \phi
  //

  double t = atof(argv[2]);
  double dphi = 2*M_PI / 60.; 
  double xl[6] = {0,     (M_PI/3.) - 0.002967059,  0-0.002967059, 0,     M_PI/3.-0.002967059, dphi*t - 0.002967059 }; 
  double xu[6] = {40000, (M_PI/3.) + 0.002967059,  0+0.002967059, 40000, M_PI/3.+0.002967059, dphi*t + 0.002967059 };

  const gsl_rng_type *T;
  gsl_rng *r;

  //gsl_monte_function G = { &gaussian, 3 , 0 };
  gsl_monte_function G = {&C00_t, 6 , NULL};

  size_t calls = atoi(argv[1]);

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (6);
  gsl_monte_vegas_integrate (&G, xl, xu, 6, 10000, r, s, &res, &err);

  do{
    gsl_monte_vegas_integrate (&G, xl, xu, 6, calls/5, r, s, &res, &err);
    //std::cout << (gsl_monte_vegas_chisq(s) -1.0) << std::endl;
  }  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
  	
  printf("%f\t%f\n",res,err);
  gsl_monte_vegas_free (s);
  

  gsl_rng_free (r);

  return 0;
}

