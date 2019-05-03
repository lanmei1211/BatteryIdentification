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


#define N      400   /* number of data points to fit */
#define TMAX   (40.0) /* time variable in [0,TMAX] */
#define MAXW 10


double w_arr[N/2];
double Z_r_arr[N/2];
double Z_i_arr[N/2];

struct data
{
    size_t n;
    double * t;
    double * y;
};

double chisq, chisq0;
    int status, info;
    size_t i;

    //Servono per il setup dell'algoritmo
    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol;
    gsl_multifit_nlinear_type *T;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params;
    const size_t n = N;
    const size_t p = 9;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar;
    double t[N], y[N], weights[N];
    struct data d;
    double x_init[9] = { 5e-6,0.1,0.77,0.6,0.2,0.60,0.7,1,0.18};
                                            /* starting values */
    //double x_init[9] = { 4e-06,0.2,0.88,0.4,0.1,0.70,0.7,1,0.18};
    gsl_vector_view x;
    gsl_vector_view wts;
    gsl_rng * r;
//Bound da impostare con min,max
//double L_bounds[2] = {3e-6,8e-6};
//double Rm_bounds[2]={0,0.3};
//double Q1_bounds[2]={0,0.80};
//double a1_bounds[2]={0.5,1};
//double Rp1_bounds[2]={0,1};
//double Q2_bounds[2]={0,0.70};
//double a2_bounds[2]={0.5,1};
//double Rp2_bounds[2] = {0,1};
//double Aw_bounds[2] = {0,0.50};



int batteria(const gsl_vector * x, void *data,
         gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *t = ((struct data *)data)->t;
    double *y = ((struct data *)data)->y;

    double L = gsl_vector_get (x, 0);
    double Rm = gsl_vector_get (x, 1);
    double Q1 = gsl_vector_get (x, 2);
    double a1 = gsl_vector_get (x, 3);
    double Rp1 = gsl_vector_get (x, 4);
    double Q2 = gsl_vector_get (x, 5);
    double a2 = gsl_vector_get (x, 6);
    double Rp2 = gsl_vector_get (x, 7);
    double Aw = gsl_vector_get (x, 8);

    //Formula per applicare i bound
//    L=L_bounds[0]+(L_bounds[1]-L_bounds[0])/(1+exp(-L));
//    Rm=Rm_bounds[0]+(Rm_bounds[1]-Rm_bounds[0])/(1+exp(-Rm));
//    Q1=Q1_bounds[0]+(Q1_bounds[1]-Q1_bounds[0])/(1+exp(-Q1));                    //e^-oo fa 0, e^oo fa infinito
//    a1=a1_bounds[0]+(a1_bounds[1]-a1_bounds[0])/(1+exp(-a1));
//    Rp1=Rp1_bounds[0]+(Rp1_bounds[1]-Rp1_bounds[0])/(1+exp(-Rp1));
//    Q2=Q2_bounds[0]+(Q2_bounds[1]-Q2_bounds[0])/(1+exp(-Q2));
//    a2=a2_bounds[0]+(a2_bounds[1]-a2_bounds[0])/(1+exp(-a2));
//    Rp2=Rp2_bounds[0]+(Rp2_bounds[1]-Rp2_bounds[0])/(1+exp(-Rp2));
//    Aw=Aw_bounds[0]+(Aw_bounds[1]-Aw_bounds[0])/(1+exp(-Aw));


    size_t i;
    //int temp = 0;
    for (i = 0; i < n; i++)
    {
        double w = t[i];
        double  Yi;
        //vars: Rm, Rp, a, Q
        //double complex battery_impedance = Rm + 1/(Q*cpow(w*I,a)+1/(Rp+Zw));
        //double complex battery_impedance = L*w*I + Rm + 1/(Q1*cpow(w*I,a1)+1/(Rp1)) +1/(Q2*cpow(w*I,a2)+1/(Rp2)) + (Aw/sqrt(w))*(1-I) ;
        double complex battery_impedance = L*I*w + Rm + 1/( 1/(1/(Q1*cpow(I*w,a1))) + 1/Rp1 ) +  1/( 1/(1/(Q2*cpow(I*w,a2))) + 1/Rp2 ) + (Aw*(1 - I))/sqrt(w);

        // Complex formato di gsl , parte im e real z = x+iy
        // x è double y è double
        // z complex
        if(i<N/2){
                //printf("REAL");
        //double re = R0 + R1/(1+pow(C1,2)*pow(w,2)*(pow(R1,2)))+ R2/(1+pow(C2,2)*pow(w,2)*(pow(R2,2)));
            double re = creal(battery_impedance);
            Yi = re;
        }else{
            //printf("IMG");
        //double im = -(pow(R1,2)*w*C1)/(1+pow(C1,2)*pow(w,2)*pow(R1,2)) - (pow(R2,2)*w*C2)/(1+pow(C2,2)*pow(w,2)*(pow(R2,2)));
            //temp = 1;
            double im = cimag(battery_impedance);
            Yi = im;
        }
        double weight = 1;
        if(a1>1 || a1<0.4){
            weight = MAXW;
        }
         if(a2>1 || a2<0.4){
            weight = MAXW;
        }

        if(Rm<0 || Rm>0.3){
            weight = MAXW;
        }

        if(Rp1<0 || Rp1 >1){
            weight = MAXW;
        }
        if(Rp2<0 || Rp2 > 2){
            weight = MAXW;
        }

        if(L < 1e-6 || L>8e-6){
            weight = MAXW;
        }

        if(Q1<0 || Q1 > 0.8){
             weight = MAXW;
        }
        if(Q2<0 || Q2 > 0.8){
             weight = MAXW;
        }
        if(Aw<0.01 || Aw > 0.5){
             weight = MAXW;
        }
        gsl_vector_set (f, i,weight*(Yi - y[i]));
}
    return GSL_SUCCESS;
}
void leggiFile(double *dati, char* nomefile){
    //Funzione che legge da file binario dei dati e li inserisce nel vettore passato come argomento
    FILE *fd;
    //double pulsazioni[100];
    int res;

    //fd=fopen("wZdot01.dat", "rb");
    fd=fopen(nomefile, "rb");
    if( fd ==NULL ) {
        perror("Errore in apertura del file");
        exit(1);
    }

    res = fread(dati, sizeof(double), 200, fd);
    fclose(fd);
}

void write_files(double *dati,char* nomefile){

    FILE *fd;
    int res;
    int k ;
    fd = fopen(nomefile,"w");
    if( fd ==NULL ) {
        perror("Errore in apertura del file");
        exit(1);
    }


    for(k = 0;k<9;k++){
            fprintf(fd,"%.15e ",dati[k]);
        }

    //res = fwrite(dati,sizeof(double),9,fd);
    //fprintf(fd,dati[0]);
    fclose(fd);

}


void read_files(){

    int contatore;
    leggiFile(w_arr,"w_bat1.bin");
    //leggiFile(w_arr,"w_modello6.bin");
    for(contatore = 0; contatore < 200; contatore++){
        //precisione matlab
        printf("%d: %.15e \n",contatore,w_arr[contatore]);
    }
    leggiFile(Z_r_arr,"real_bat1.bin");
    leggiFile(Z_i_arr,"imag_bat1.bin");
    //leggiFile(Z_r_arr,"real_modello6.bin");
    //leggiFile(Z_i_arr,"imag_modello6.bin");
    }

void callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    double rcond;

    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stderr, "iter %2zu: L = %.e, Rm = %.e, Q1 = %.e, a1 = %.e, Rp1 = %.e, Q2 = %.e,a2 = %.e, Rp2 = %.e, Aw = %.e, cond(J) = %8.4f, |f(x)| = %.4f\n",
            iter,
            gsl_vector_get(x, 0),
            gsl_vector_get(x, 1),
            gsl_vector_get(x, 2),
            gsl_vector_get(x, 3),
            gsl_vector_get(x, 4),
            gsl_vector_get(x, 5),
            gsl_vector_get(x, 6),
            gsl_vector_get(x, 7),
            gsl_vector_get(x, 8),
            1.0 / rcond,
            gsl_blas_dnrm2(f));
}


void identification(){





    initialize_params();
    acquire_data();
    solve_system();
    print_results();

}
void initialize_params(){
T = gsl_multifit_nlinear_trust;
    fdf_params = gsl_multifit_nlinear_default_parameters();
    covar = gsl_matrix_alloc (p, p);
    gsl_rng_env_setup();
    d.n=n;
    d.t = t;
    d.y=y;
     x = gsl_vector_view_array (x_init, p);
     wts = gsl_vector_view_array(weights, n);
    r = gsl_rng_alloc(gsl_rng_default);
    read_files();

    /* define the function to be minimized */
    fdf.f = batteria;      //EQ DA CAMBIARE
    fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = n;
    fdf.p = p;
    fdf.params = &d;

}
void acquire_data(){

    /* this is the data to be fitted */
    for (i = 0; i < n; i++)
    {

        if(i<200){
           t[i] = w_arr[i];
           y[i] = Z_r_arr[i];


        } else {
           t[i] = w_arr[i-200];
           y[i] = Z_i_arr[i-200];

        }

        //double si = 0.1 * y[i];
        //weights[i] = 1.0 / (si * si);

        //DA GUARDARE per il nuovo controllo
        weights[i] = 1.0;
        printf ("data: %g %g\n", t[i], y[i]);
    };

}

void solve_system(){

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system with a maximum of 100 iterations */
    status = gsl_multifit_nlinear_driver(5000, xtol, gtol, ftol,
                                         callback, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);


}

void print_results(){

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

    {
        double dof = n - p;
        double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        //double L=L_bounds[0]+(L_bounds[1]-L_bounds[0])/(1+exp(-FIT(0)));
        double L = FIT(0);
        //double Rm=Rm_bounds[0]+(Rm_bounds[1]-Rm_bounds[0])/(1+exp(-FIT(1)));
        double Rm = FIT(1);
        //double Q1=Q1_bounds[0]+(Q1_bounds[1]-Q1_bounds[0])/(1+exp(-FIT(2)));                    //e^-oo fa 0, e^oo fa infinito
        double Q1 = FIT(2);
        //double a1=a1_bounds[0]+(a1_bounds[1]-a1_bounds[0])/(1+exp(-FIT(3)));
        double a1 = FIT(3);
        //double Rp1=Rp1_bounds[0]+(Rp1_bounds[1]-Rp1_bounds[0])/(1+exp(-FIT(4)));
        double Rp1 = FIT(4);
        //double Q2=Q2_bounds[0]+(Q2_bounds[1]-Q2_bounds[0])/(1+exp(-FIT(5)));
        double Q2 = FIT(5);
        //double a2=a2_bounds[0]+(a2_bounds[1]-a2_bounds[0])/(1+exp(-FIT(6)));
        double a2 = FIT(6);
        //double Rp2=Rp2_bounds[0]+(Rp2_bounds[1]-Rp2_bounds[0])/(1+exp(-FIT(7)));
        double Rp2 = FIT(7);
        //double Aw=Aw_bounds[0]+(Aw_bounds[1]-Aw_bounds[0])/(1+exp(-FIT(8)));
        double Aw = FIT(8);

        double dati[9] = {L,Rm,Q1,a1,Rp1,Q2,a2,Rp2,Aw};
        write_files(dati,"res_bat1.dat");

        fprintf (stderr, "L      = %.15e +/- %.15e\n", L, c*ERR(0));
        fprintf (stderr, "Rm     = %.15e +/- %.15e\n", Rm, c*ERR(1));
        fprintf (stderr, "Q1      = %.15e +/- %.15e\n", Q1, c*ERR(2));
        fprintf (stderr, "a1      = %.15e +/- %.15e\n", a1, c*ERR(3));
        fprintf (stderr, "Rp1     = %.15e +/- %.15e\n", Rp1, c*ERR(4));
        fprintf (stderr, "Q2     = %.15e +/- %.15e\n", Q2, c*ERR(5));
        fprintf (stderr, "a2      = %.15e +/- %.15e\n", a2, c*ERR(6));
        fprintf (stderr, "Rp2      = %.15e +/- %.15e\n", Rp2, c*ERR(7));
        fprintf (stderr, "Aw      = %.15e +/- %.15e\n", Aw, c*ERR(8));
    }

    fprintf (stderr, "status = %s\n", gsl_strerror (status));


}

void close(){
    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);
}

int main (void)
{
    identification();
    return 0;
}
