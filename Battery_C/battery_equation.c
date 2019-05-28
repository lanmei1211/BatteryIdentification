#include "battery_equation.h"
#include "global.h"

int batteria(const gsl_vector * x, void *data, gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *t = ((struct data *)data)->t;
    double *y = ((struct data *)data)->y;

    size_t i;

    for (i = 0; i < n; i++)
    {
        double w = t[i];
        double  Yi;


        double complex battery_impedance;
        if(type == 1){
            battery_impedance = R_RCPE(w,P(0),P(1),P(2),P(3));
        } else if (type == 2){
            battery_impedance = R_RC_RC(w,P(0),P(1),P(2),P(3),P(4));
        } else if (type == 3){
            battery_impedance = CPE_W(w,P(0),P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8));
        } else {
            battery_impedance = fouquet(w,P(0),P(1),P(2),P(3),P(4),P(5));
        }

        if(i<N/2){
            double re = creal(battery_impedance);
            Yi = re;
        }else{
            double im = cimag(battery_impedance);
            Yi = im;
        }

        double weight = check_limit(x);
        gsl_vector_set (f, i,weight*(Yi - y[i]));
    }
    return GSL_SUCCESS;
}

double check_limit(gsl_vector * x){
    for(int i=0;i<p;i++){
        double value = P(i);
        double ub = upper_bounds[i];
        double lb = lower_bounds[i];
        if (value < lb || value > ub){
            return MAXW;
        }
    }
    return 1;
}

double complex R_RCPE(double w,double Rs, double Q1, double n1, double Rp1){
    return Rs+1/( 1/(1/(Q1*cpow(I*w,n1))) + 1/Rp1);
}

double complex R_RC_RC(double w, double Rs,double C1,double Rp1,double C2,double Rp2){
    return Rs + 1/( 1/(1/(C1*(I*w))) + 1/Rp1 ) +  1/( 1/(1/(C2*(I*w))) + 1/Rp2 );
}

double complex CPE_W(double w, double L, double Rs, double Q1, double n1, double Rp1, double Q2, double n2, double Rp2, double Aw){
    return L*I*w + Rs + 1/( 1/(1/(Q1*cpow(I*w,n1))) + 1/Rp1 ) +  1/( 1/(1/(Q2*cpow(I*w,n2))) + 1/Rp2 ) + (Aw*(1 - I))/sqrt(w);

}

double complex fouquet(double w, double a, double b, double c, double d, double e, double f){
    return 0;
}
