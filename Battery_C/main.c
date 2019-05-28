#include "main.h"
#include "classifier.h"
#include "battery_equation.h"
#include "global.h"
#include "filter.h"

double w_arr[N/2];
double Z_r_arr[N/2];
double Z_i_arr[N/2];
double chisq, chisq0;
int status, info;
int i;
gsl_multifit_nlinear_type *T;
gsl_multifit_nlinear_workspace *w;
gsl_multifit_nlinear_fdf fdf;
gsl_multifit_nlinear_parameters fdf_params;


gsl_vector *f;
gsl_matrix *J;
gsl_matrix *covar;
double t[N], y[N], weights[N];
Data d;

gsl_vector_view x;
gsl_rng * r;

int type;
size_t p;
double x_init[9];
double lower_bounds[9];
double upper_bounds[9];

int fuel_cell = 0;


void leggiFile(double *dati, char* nomefile){
    //Funzione che legge da file binario dei dati e li inserisce nel vettore passato come argomento
    FILE *fd;

    int res;

    fd=fopen(nomefile, "rb");
    if( fd ==NULL ) {
        perror("Errore in apertura del file");
        exit(1);
    }

    res = fread(dati, sizeof(double), N/2, fd);
    fclose(fd);
}

void write_files(double *dati,char* nomefile,int n){

    FILE *fd;
    int res;
    int k ;
    fd = fopen(nomefile,"w");
    if( fd ==NULL ) {
        perror("Errore in apertura del file");
        exit(1);
    }


    for(k = 0;k<n;k++){
            //fprintf(fd,"%.15e ",dati[k]);
            fwrite(&dati[k],sizeof(double),1,fd);
        }


    fclose(fd);

}


void read_files(){

    int contatore;
    leggiFile(w_arr,"w_bat1.bin");
    for(contatore = 0; contatore < N/2; contatore++){
        //precisione matlab
        printf("%d: %.15e \n",contatore,w_arr[contatore]);
    }
    leggiFile(Z_r_arr,"real_bat1.bin");
    leggiFile(Z_i_arr,"imag_bat1.bin");
}

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    double rcond;

    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&rcond, w);

    /*
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
    */
}


void identification(){

    printf("Acquire data\n\n");
    acquire_data();
    printf("Classify\n\n");
    classify();
    printf("Initialize Params\n\n");
    initialize_params();
    printf("Solve System\n\n");
    solve_system();
    printf("Print results\n\n");
    print_results();
    printf("Close\n\n");
    close();

}
void initialize_params(){
    T = gsl_multifit_nlinear_trust;
    fdf_params = gsl_multifit_nlinear_default_parameters();
    covar = gsl_matrix_alloc (p, p);
    gsl_rng_env_setup();
    d.n=N;
    d.t = t;
    d.y=y;
    x = gsl_vector_view_array (x_init, p);
    r = gsl_rng_alloc(gsl_rng_default);


    /* define the function to be minimized */
    fdf.f = batteria;
    fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = N;
    fdf.p = p;
    fdf.params = &d;

}
void acquire_data(){
    read_files();
    /* this is the data to be fitted */
    for (i = 0; i < N; i++)
    {

        if(i<N/2){
           t[i] = w_arr[i];
           y[i] = Z_r_arr[i];


        } else {
           t[i] = w_arr[i-N/2];
           y[i] = Z_i_arr[i-N/2];

        }


        printf ("data: %g %g\n", t[i], y[i]);
    };

}

void classify(){
    if(fuel_cell == 0){
        double filtered_data[N/2];
        filter(Z_i_arr,filtered_data,N/2);
        type = classificatore(filtered_data+OFFSET,N/2-OFFSET);
    } else {
        type = 4;
    }
    printf("Type is: %d\n", type);
    if (type == 1){
        p=4;

        x_init[0] = 0.1;
        x_init[1] = 0.77;
        x_init[2] = 0.6;
        x_init[3] = 0.2;

        lower_bounds[0] = 0.05;
        lower_bounds[1] = 0.74;
        lower_bounds[2] = 0.48;
        lower_bounds[3] = 0.08;

        upper_bounds[0] = 0.15;
        upper_bounds[1] = 0.80;
        upper_bounds[2] = 0.72;
        upper_bounds[3] = 0.24;

    } else if (type == 2){
        p=5;

        x_init[0] = 0.15;
        x_init[1] = 1;
        x_init[2] = 0.2;
        x_init[3] = 5;
        x_init[4] = 1;

        lower_bounds[0] = 0.1;
        lower_bounds[1] = 1;
        lower_bounds[2] = 0.2;
        lower_bounds[3] = 5;
        lower_bounds[4] = 1;


        upper_bounds[0] = 0.2;
        upper_bounds[1] = 1.2;
        upper_bounds[2] = 0.24;
        upper_bounds[3] = 5.1;
        upper_bounds[4] = 1.1;


    } else if (type == 3){
        p=9;

        x_init[0] = 3e-6;
        x_init[1] = 0.2;
        x_init[2] = 0.04;
        x_init[3] = 0.6;
        x_init[4] = 0.4;
        x_init[5] = 0.6;
        x_init[6] = 0.8;
        x_init[7] = 0.7;
        x_init[8] = 0.22;

        lower_bounds[0] = 2e-6;
        lower_bounds[1] = 0.0;
        lower_bounds[2] = 0.0;
        lower_bounds[3] = 0.5;
        lower_bounds[4] = 0.0;
        lower_bounds[5] = 0.0;
        lower_bounds[6] = 0.5;
        lower_bounds[7] = 0.0;
        lower_bounds[8] = 0.0;

        upper_bounds[0] = 8e-6;
        upper_bounds[1] = 0.3;
        upper_bounds[2] = 0.8;
        upper_bounds[3] = 1.0;
        upper_bounds[4] = 1.0;
        upper_bounds[5] = 0.7;
        upper_bounds[6] = 1.0;
        upper_bounds[7] = 1.0;
        upper_bounds[8] = 0.5;
    } else {
        p=9;

        x_init[0] = 3e-6;
        x_init[1] = 0.2;
        x_init[2] = 0.04;
        x_init[3] = 0.6;
        x_init[4] = 0.4;
        x_init[5] = 0.6;

        lower_bounds[0] = 2e-6;
        lower_bounds[1] = 0.0;
        lower_bounds[2] = 0.0;
        lower_bounds[3] = 0.5;
        lower_bounds[4] = 0.0;
        lower_bounds[5] = 0.0;

        upper_bounds[0] = 8e-6;
        upper_bounds[1] = 0.3;
        upper_bounds[2] = 0.8;
        upper_bounds[3] = 1.0;
        upper_bounds[4] = 1.0;
        upper_bounds[5] = 0.7;

    }

}


void solve_system(){

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, N, p);

    /* initialize solver with starting point */
    gsl_multifit_nlinear_init (&x.vector, &fdf, w);

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
        double dof = N - p;
        double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        if (type==1){
            double dati[4] = {FIT(0),FIT(1),FIT(2),FIT(3)};
            write_files(dati,"res_bat1.dat",p);
        }else if(type==2){
            double dati[5] = {FIT(0),FIT(1),FIT(2),FIT(3),FIT(4)};
            write_files(dati,"res_bat1.dat",p);
        } else if(type == 3){
            double dati[9] = {FIT(0),FIT(1),FIT(2),FIT(3),FIT(4),FIT(5),FIT(6),FIT(7),FIT(8)};
            write_files(dati,"res_bat1.dat",p);
        } else {
            double dati[6] = {FIT(0),FIT(1),FIT(2),FIT(3),FIT(4),FIT(5)};
            write_files(dati,"res_bat1.dat",p);
        }



        fprintf (stderr, "Param 1      = %.15e \n", FIT(0));
        fprintf (stderr, "Param 2     = %.15e \n", FIT(1));
        fprintf (stderr, "Param 3      = %.15e \n", FIT(2));
        fprintf (stderr, "Param 4      = %.15e \n", FIT(3));
        if(type!=1){
            fprintf (stderr, "Param 5     = %.15e \n", FIT(4));
            if(type==4){
                fprintf (stderr, "Param 6     = %.15e \n", FIT(5));
            }
            if(type==3){
                fprintf (stderr, "Param 6     = %.15e \n", FIT(5));
                fprintf (stderr, "Param 7      = %.15e \n", FIT(6));
                fprintf (stderr, "Param 8      = %.15e \n", FIT(7));
                fprintf (stderr, "Param 9      = %.15e  \n", FIT(8));
            }
        }
    }

    fprintf (stderr, "status = %s\n", gsl_strerror (status));


}

void close(){
    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);
}

int main (int argc, char **argv)
{
    if(argc ==2){
        fuel_cell=1;
    }
    //test_write();
    identification();
    return 0;
}
