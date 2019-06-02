/**
* Main.c
*/

#include "main.h"
#include "classifier.h"
#include "battery_equation.h"
#include "global.h"
#include "filter.h"

/**
* numero di punti rilevato
* rappresenta il numero di frequenze disponibili
* viene inizializzato dinamicamente analizzando il file dei dati
*/
int N;

double *w_arr;
double *Z_r_arr;
double *Z_i_arr;
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
double *t, *y, *weights;
Data d;

gsl_vector_view x;
gsl_rng * r;

/**
* tipo di modello da utilizzare
*/
int type;
/**
* numero di parametri
*/
size_t p;
double x_init[9];
double lower_bounds[9];
double upper_bounds[9];

/**
* indica se si sta analizzando una fuel cell o una batteria a litio
*/
int fuel_cell = 0;

/**
* aumenta dimensione dell'area di memoria dedicata a un puntatore
* @param ptr puntatore
* @param elements valori da salvare
*/
int increase_size(double** ptr,int elements){
    *ptr =realloc(*ptr, elements*sizeof(double));
}

/**
* leggi i dati da un file
* @param dati variabile in cui salvare i dati
* @param nomefile nome del file
*/
void read_single_file(double **dati, char* nomefile){
    FILE *fd;

    fd=fopen(nomefile, "rb");
    if( fd ==NULL ) {
        perror("Errore in apertura del file");
        exit(1);
    }

    int elements =0;
    *dati = malloc(sizeof(double)*(elements+1));
    while(fread((*dati)+elements, sizeof(double), 1, fd) > 0){
        elements++;
        increase_size(dati, elements+1);


    }
    N=elements*2;
    fclose(fd);
}

/**
* scrivi il risultato dell'identificazione su file binario
* @param dati valori da scrivere
* @param nomefile nome del file
* @param n numero di parametri
*/
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
            fwrite(&dati[k],sizeof(double),1,fd);
        }


    fclose(fd);

}

/**
* leggi i file contenenti i dati
* @see read_single_file()
*/
void read_files(){

    read_single_file(&w_arr,"w_bat1.bin");
    read_single_file(&Z_r_arr,"real_bat1.bin");
    read_single_file(&Z_i_arr,"imag_bat1.bin");
}

/**
* funzione richiamata dopo una iterazione del curve fitter
*/
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

/**
* identifica un insieme di dati di una batteria
* @see acquire_data()
* @see classify()
* @see initialize_params()
* @see solve_system()
* @see print_result()
* @see close()
*/
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

/**
* inizializza l'identificatore
*/
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
/**
* acquisisci i dati da file e inseriscili negli array utilizzati da gsl
* @see read_files()
*/
void acquire_data(){
    read_files();
    t = malloc(sizeof(double)*N);
    y = malloc(sizeof(double)*N);
    weights = malloc(sizeof(double)*N);
    /* this is the data to be fitted */
    printf("N vale %d\n",N);
    for (i = 0; i < N; i++)
    {

        if(i<N/2){
           t[i] = w_arr[i];
           y[i] = Z_r_arr[i];


        } else {
           t[i] = w_arr[i-N/2];
           y[i] = Z_i_arr[i-N/2];

        }


        //printf ("data: %g %g\n", t[i], y[i]);
    };

}
/**
* classifica il tipo di curva, se i dati sono per una fuel cell, assegna direttamente i parametri,
* altrimenti, filtra i dati e chiama il classificatore
* @see filter()
* @see classificatore()
*/
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

        lower_bounds[0] = 0.02;
        lower_bounds[1] = 0.72;
        lower_bounds[2] = 0.48;
        lower_bounds[3] = 0.02;

        upper_bounds[0] = 0.22;
        upper_bounds[1] = 0.82;
        upper_bounds[2] = 0.72;
        upper_bounds[3] = 0.3;

    } else if (type == 2){
        p=5;

        x_init[0] = 0.15;
        x_init[1] = 1;
        x_init[2] = 0.2;
        x_init[3] = 5;
        x_init[4] = 1;

        lower_bounds[0] = 0.07;
        lower_bounds[1] = 0.8;
        lower_bounds[2] = 0.15;
        lower_bounds[3] = 4.7;
        lower_bounds[4] = 0.9;


        upper_bounds[0] = 0.25;
        upper_bounds[1] = 1.2;
        upper_bounds[2] = 0.24;
        upper_bounds[3] = 5.3;
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
        p=6;

        x_init[0] = 0.005;
        x_init[1] = 1;
        x_init[2] = 0.7;
        x_init[3] = 0.01;
        x_init[4] = 0.05;
        x_init[5] = 0.01;

        lower_bounds[0] = 0.0;
        lower_bounds[1] = 0.5;
        lower_bounds[2] = 0.5;
        lower_bounds[3] = 0.0;
        lower_bounds[4] = 0.0;
        lower_bounds[5] = 0.0;

        upper_bounds[0] = 1.0;
        upper_bounds[1] = 2.0;
        upper_bounds[2] = 1.0;
        upper_bounds[3] = 1.0;
        upper_bounds[4] = 1.0;
        upper_bounds[5] = 0.5;

    }

}

/**
* risolvi il sistema di identificazione parametrica
*/
void solve_system(){

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, N, p);

    /* initialize solver with starting point */
    gsl_multifit_nlinear_init (&x.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);


    /* solve the system with a maximum of 5000 iterations */
    status = gsl_multifit_nlinear_driver(5000, xtol, gtol, ftol,
                                         callback, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);


}

/**
* stampa risultati dell'identificazione, sia a video che su file 'res_bat1.dat'
* @see write_files()
*/
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
        char file_name[] = "res_bat1.dat";
        if (type==1){
            double dati[4] = {FIT(0),FIT(1),FIT(2),FIT(3)};
            write_files(dati,file_name,p);
        }else if(type==2){
            double dati[5] = {FIT(0),FIT(1),FIT(2),FIT(3),FIT(4)};
            write_files(dati,file_name,p);
        } else if(type == 3){
            double dati[9] = {FIT(0),FIT(1),FIT(2),FIT(3),FIT(4),FIT(5),FIT(6),FIT(7),FIT(8)};
            write_files(dati,file_name,p);
        } else {
            double dati[6] = {FIT(0),FIT(1),FIT(2),FIT(3),FIT(4),FIT(5)};
            write_files(dati,file_name,p);
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

/**
* libera lo spazio allocato per l'identificazione
*/
void close(){
    free(w_arr);
    free(Z_r_arr);
    free(Z_i_arr);
    free(t);
    free(y);
    free(weights);
    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);
}
/**
* main del programma di identificazione parametri
* deve avere due argomenti, di cui il secondo deve essere 1 per fuel cell e 0 per batterie a litio
* @param argc numero di argomenti
* @param argv argomenti
* @see identification();
* @return 0
*/

int main (int argc, char **argv)
{
    if(argc == 2){

        if(strcmp(argv[1],"1") == 0){
            fuel_cell = 1;
        } else if(strcmp(argv[1],"0") == 0){
            fuel_cell = 0;
        } else {
            printf("Errore, selezionare il tipo di batteria:\n1 per Fuel cell, 0 per batteria al litio\n");
            return 0;
        }
        identification();
    } else {
        printf("Errore, selezionare il tipo di batteria:\n1 per Fuel cell, 0 per batteria al litio\n");
    }
    return 0;
}
