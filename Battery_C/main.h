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
#include <complex.h>


#define N 400
#define MAXW 10
const size_t p = 9;
const double xtol = 1e-8;
const double gtol = 1e-8;
const double ftol = 0.0;
const double x_init[9] = { 5e-6,0.1,0.77,0.6,0.2,0.60,0.7,1,0.18};


typedef struct data
{
    size_t n;
    double * t;
    double * y;
} Data;


int batteria(const gsl_vector *, void *,gsl_vector *);
void leggiFile(double*, char* );
void write_files(double *,char* );
void read_files();
void callback(const size_t , void *,const gsl_multifit_nlinear_workspace *);
void identification();
void initialize_params();
void acquire_data();
void solve_system();
void print_results();
void close();

#endif // MAIN_H_INCLUDED
