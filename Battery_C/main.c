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

struct timespec start, end;

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

char base_path[] = "data/modelloFuelZdot1/";

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


uint64_t elapsed_time(){
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    return delta_us;
}


/**
* aumenta dimensione dell'area di memoria dedicata a un puntatore
* @param ptr puntatore
* @param elements valori da salvare
*/
int increase_size(double** ptr,int elements){
    *ptr =realloc(*ptr, elements*sizeof(double));
    return 0;
}

/**
* leggi i dati da un file
* @param dati variabile in cui salvare i dati
* @param nomefile nome del file
*/
void read_single_file(double **dati, char* nomefile){
    char *path = malloc(strlen(base_path) + strlen(nomefile)+1);
    strcpy(path,base_path);
    strcat(path,nomefile);
    FILE *fd;

    fd=fopen(path, "rb");
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
    free(path);
}

/**
* scrivi il risultato dell'identificazione su file binario
* @param dati valori da scrivere
* @param nomefile nome del file
* @param n numero di parametri
*/
void write_files(double *dati,char* nomefile,int n){
    char *path = malloc(strlen(base_path) + strlen(nomefile)+1);
    strcpy(path,base_path);
    strcat(path,nomefile);

    FILE *fd;
    int k ;
    fd = fopen(path,"w");
    if( fd ==NULL ) {
        perror("Errore in apertura del file");
        exit(1);
    }


    for(k = 0;k<n;k++){
            fwrite(&dati[k],sizeof(double),1,fd);
        }


    fclose(fd);
    free(path);

}

/**
* leggi i file contenenti i dati
* @see read_single_file()
*/
void read_files(){

    read_single_file(&w_arr,"w.bin");
    read_single_file(&Z_r_arr,"real.bin");
    read_single_file(&Z_i_arr,"imag.bin");
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

    clock_gettime(CLOCK_MONOTONIC, &end);
    uint64_t elapsed = elapsed_time();
    printf("Iteration number: %2zu, Elapsed time: %lu\n",iter, elapsed);
    clock_gettime(CLOCK_MONOTONIC, &start);

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


        printf ("data %d: %g %g\n",i, t[i], y[i]);
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
        type = FUEL_CELL;
    }
    printf("Type is: %d\n", type);
    if (type == ONE_CURVE){
        p=4;

        x_init[0] = 0.085;
        x_init[1] = 0.72;
        x_init[2] = 0.5;
        x_init[3] = 0.13;

        lower_bounds[0] = 0.07;
        lower_bounds[1] = 0.60;
        lower_bounds[2] = 0.45;
        lower_bounds[3] = 0.02;

        upper_bounds[0] = 0.15;
        upper_bounds[1] = 0.90;
        upper_bounds[2] = 0.80;
        upper_bounds[3] = 0.3;

    } else if (type == TWO_CURVES){
        p=5;

        x_init[0] = 0.085;
        x_init[1] = 0.72;
        x_init[2] = 0.5;
        x_init[3] = 4.20;
        x_init[4] = 1.1;

        lower_bounds[0] = 0.07;
        lower_bounds[1] = 0.60;
        lower_bounds[2] = 0.1;
        lower_bounds[3] = 3.80;
        lower_bounds[4] = 0.9;


        upper_bounds[0] = 0.15;
        upper_bounds[1] = 1.20;
        upper_bounds[2] = 0.70;
        upper_bounds[3] = 6.20;
        upper_bounds[4] = 1.2;


    } else if (type == TAIL){
        p=9;

        x_init[0] = 4e-6;
        x_init[1] = 0.05;
        x_init[2] = 0.025;
        x_init[3] = 0.68;
        x_init[4] = 0.43;
        x_init[5] = 0.3;
        x_init[6] = 0.8;
        x_init[7] = 0.6;
        x_init[8] = 0.22;

        lower_bounds[0] = 3e-6;
        lower_bounds[1] = 0.02;
        lower_bounds[2] = 0.0;
        lower_bounds[3] = 0.5;
        lower_bounds[4] = 0.2;
        lower_bounds[5] = 0.1;
        lower_bounds[6] = 0.5;
        lower_bounds[7] = 0.3;
        lower_bounds[8] = 0.1;

        upper_bounds[0] = 6e-6;
        upper_bounds[1] = 0.07;
        upper_bounds[2] = 0.08;
        upper_bounds[3] = 0.8;
        upper_bounds[4] = 0.7;
        upper_bounds[5] = 0.7;
        upper_bounds[6] = 1.0;
        upper_bounds[7] = 0.9;
        upper_bounds[8] = 0.5;
    } else {
        p=6;
        //modelloFuelZdot1

        x_init[0] = 0.006;
        x_init[1] = 1.9;
        x_init[2] = 0.85;
        x_init[3] = 0.002;
        x_init[4] = 0.02;
        x_init[5] = 0.6;

        lower_bounds[0] = 0.001;
        lower_bounds[1] = 1.8;
        lower_bounds[2] = 0.7;
        lower_bounds[3] = 0.002;
        lower_bounds[4] = 0.01;
        lower_bounds[5] = 0.2;

        upper_bounds[0] = 0.01;
        upper_bounds[1] = 2.3;
        upper_bounds[2] = 1.2;
        upper_bounds[3] = 0.01;
        upper_bounds[4] = 0.3;
        upper_bounds[5] = 0.7;


        //Zcella15
        /*
        x_init[0] = 0.04;
        x_init[1] = 0.25;
        x_init[2] = 0.65;
        x_init[3] = 0.08;
        x_init[4] = 0.06;
        x_init[5] = 0.2;

        lower_bounds[0] = 0.02;
        lower_bounds[1] = 0.1;
        lower_bounds[2] = 0.4;
        lower_bounds[3] = 0.02;
        lower_bounds[4] = 0.005;
        lower_bounds[5] = 0.07;

        upper_bounds[0] = 0.07;
        upper_bounds[1] = 0.5;
        upper_bounds[2] = 0.9;
        upper_bounds[3] = 0.4;
        upper_bounds[4] = 0.09;
        upper_bounds[5] = 0.4;
        */


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

    clock_gettime(CLOCK_MONOTONIC, &start);
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

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);
        char file_name[] = "res_pc.dat";
        if (type== ONE_CURVE){
            double dati[4] = {FIT(0),FIT(1),FIT(2),FIT(3)};
            write_files(dati,file_name,p);
        }else if(type== TWO_CURVES){
            double dati[5] = {FIT(0),FIT(1),FIT(2),FIT(3),FIT(4)};
            write_files(dati,file_name,p);
        } else if(type == TAIL){
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
        if(type!=ONE_CURVE){
            fprintf (stderr, "Param 5     = %.15e \n", FIT(4));
            if(type==FUEL_CELL){
                fprintf (stderr, "Param 6     = %.15e \n", FIT(5));
            }
            if(type==TAIL){
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
