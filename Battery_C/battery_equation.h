#ifndef BATTERY_EQUATION_H_INCLUDED
#define BATTERY_EQUATION_H_INCLUDED

#include <gsl/gsl_vector.h>
#include <complex.h>

#define P(i) gsl_vector_get(x,i)
#define MAXW 10

int batteria(const gsl_vector *, void *,gsl_vector *);
double complex R_RCPE(double ,double , double , double , double );
double complex R_RC_RC(double , double ,double ,double ,double ,double );
double complex CPE_W(double , double , double , double , double , double , double , double , double , double );
double check_limit(gsl_vector *);


extern int type;
extern size_t p;
extern double lower_bounds[9];
extern double upper_bounds[9];
#endif // BATTERY_EQUATION_H_INCLUDED
