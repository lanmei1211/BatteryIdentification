#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_sf_lambert.h>




const double xtol = 1e-6;
const double gtol = 1e-6;
const double ftol = 1e-6;

void leggiFile(double**, char* );
void write_files(double *,char* ,int n);
void read_files();
void callback(const size_t , void *,const gsl_multifit_nlinear_workspace *);
void identification();
void initialize_params();
void acquire_data();
void solve_system();
void print_results();
void close();

#endif // MAIN_H_INCLUDED
