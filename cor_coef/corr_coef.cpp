#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include "funct.hpp"


int main (int argc, char* argv[]){

  double res, err;

  // Integration ranges xl -> lower and xu -> upper
  double xl[3] = {-10., -10., -10.}; 
  double xu[3] = { 10., 10., 10.};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = { &gaussian, 3, 0 };

  size_t calls = atoi(argv[1]);

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
  gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s, &res, &err);

  do{
    gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s, &res, &err);
  }  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
  	
  printf("%f\t%f\n",res,err);
  gsl_monte_vegas_free (s);
  

  gsl_rng_free (r);

  return 0;
}

