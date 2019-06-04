
#include "classifier.h"


double vettore[tests][size] = {
    {0, 1, 2, 3, 4 ,5},
    {5, 4, 3, 2, 1, 0},
    {0, 1, 2, 3, 2, 1},
    {3, 2, 1, 0, 1, 2},
    {0, 1, 2, 1, 0, 1},
    {0, 1, 2, 3, 4 ,3},
    {1, 0, 2, 3, 4 ,3},
    {5, 4, 3, 2, 1, 2},
    {2, 1, 2, 3, 4 ,3}
};

double vettore2[tests][size] = {
    {0, 1, 2, 3, 4 ,5},
    {5, 4, 3, 2, 1, 0},
    {0, 1, 2, 3, 2, 1},
    {3, 2, 1, 0, 1, 2},
    {0, 1, 2, 1, 0, 1},
    {0, 1, 2, 1, 2, 1},
    {5, 4, 5, 4, 3, 2},
    {1, 2, 1, 2, 1, 0}
};
int test(){


    test_is_increasing();
    test_is_decreasing();
    test_compute();
    test_classificatore();
    return 0;

}


int test_is_increasing(){
    printf("---  Test is_increasing  ---\n");
    for(int i = 0; i < tests;i++){
        printf("Test numero %d\n", i);
        for(int j=0;j<size;j++) {
            printf("%.0lf ", vettore[i][j]);
        }
        int n = sizeof(vettore[i])/sizeof(vettore[i][0]);
        int val = is_increasing(vettore[i],n);
        printf("; Out: %d\n\n", val);

    }
    return 0;
}

int test_is_decreasing(){
    printf("---  Test is_decreasing  ---\n");
    for(int i = 0; i < tests;i++){
        printf("Test numero %d\n", i);
        for(int j=0;j<size;j++) {
            printf("%.0lf ", vettore[i][j]);
        }
        int n = sizeof(vettore[i])/sizeof(vettore[i][0]);
        int val = is_decreasing(vettore[i],n);
        printf("; Out: %d\n\n", val);

    }
    return 0;
}

int test_compute(){
    printf("---  Test compute  ---\n");
    for(int i = 0; i < tests;i++){
        printf("Test numero %d\n", i);
        for(int j=0;j<size;j++) {
            printf("%.0lf ", vettore[i][j]);
        }
        int n = sizeof(vettore[i])/sizeof(vettore[i][0]);
        int val = compute(vettore[i],n);
        printf("; Out: %d\n\n", val);

    }
    return 0;
}

int test_classificatore(){
    double vettore_tipo_1[13] = {1,3,5,7,8,9,10,9,8,7,5,3,1};
    double vettore_tipo_2[31] = {1,3,5,7,9,11,12,13,14,15,16,17,18,20,22,24,26,28,30,31,32,33,34,35,36,35,34,33,31,30};
    double vettore_tipo_3[31] = {1,3,5,7,9,11,12,13,14,15,16,17,18,20,22,24,26,28,30,31,32,33,34,35,36,38,40,42,44,46,50};
    printf("---  Test classificatore  ---\n1n");
    int n = 13;

    printf("Test numero 1\n");
    for(int j=0;j<n;j++) {
        printf("%.0f ", vettore_tipo_1[j]);
    }

    int val = classificatore(vettore_tipo_1,n);
    printf("\nOut: %d\n\n", val);
    n=31;
    printf("Test numero 2\n");
    for(int j=0;j<n;j++) {
        printf("%.0f ", vettore_tipo_2[j]);
    }

    val = classificatore(vettore_tipo_2,n);
    printf("\nOut: %d\n\n", val);

    printf("Test numero 3\n");
    for(int j=0;j<n;j++) {
        printf("%.0f ", vettore_tipo_3[j]);
    }

    val = classificatore(vettore_tipo_3,n);
    printf("\nOut: %d\n\n", val);


}
/**
* Controlla se un array è crescente
* @param data array da controllare
* @param n numero di dati
* @return l'indice a cui l'array diventa decrescente, 0 se è sempre crescente
*/
int is_increasing(double *data, int n){
    for(int i = 1;i<n;i++){
        if(data[i-1]>data[i])
            return i;
    }
    return 0;
}


/**
* Controlla se un array è decrescente
* @param data array da controllare
* @param n numero di dati
* @return l'indice a cui l'array diventa crescente, 0 se è sempre decrescente
*/

int is_decreasing(double *data, int n){
    for(int i = 1;i<n;i++){
        if(data[i-1]<data[i])
            return i;
    }
    return 0;
}

/**
* analizza il numero di gobbe di un array
* @param data array da controllare
* @param n numero di dati
* @return 0 se ha una gobba, 1 se ha due gobbe
*/
int compute(double *data, int n){
    // 0 = 1 gobba, 1 = 2 gobbe

    int i_cresc = is_increasing(data,n);
    printf("i_cresc: %d\n", i_cresc);
    if(i_cresc!=0){
        int i_decr = is_decreasing(data+i_cresc,n-i_cresc);
        printf("i_decr: %d\n", i_decr);
        if(i_decr!=0){
            int i_cresc2 = is_increasing(data+i_cresc+i_decr,n-i_cresc-i_decr);
            printf("i_cresc2: %d\n", i_cresc2);
            if(i_cresc2!=0){
                return TWO_CURVES;
            }
        }
    }
    return ONE_CURVE;
}
/**
* classifica un array in base al numero di gobbe e la "coda"
* si analizza l'array fino a quando non inizia a decrescere e si calcola la pendenza discreta per ogni punto
* si controlla per i dati rimanenti se sono crescenti o decrescenti, in caso essi siano crescenti significa che
* è presente la coda, altrimenti si analizzano le pendenze calcolate precedentemente per capire il numero di gobbe
* @param data array da controllare
* @param n numero di dati
* @return 3 se ha la "coda", 2 se ha due gobbe, 1 se ha una gobba
*/
int classificatore(double *data, int n){
    int k;
    for(k=0;k<n-4;k++ ){
        printf("data class: %.04f\n", data[k]);
        if(data[k]>data[k+1]){
            k++;
            break;
        }
        data[k] = (data[k+4]-data[k])/5;
    }
    printf("Class - k: %d, n: %d\n", k,n);
    for(int j=k+1;j<n;j++){
        printf("data class: %.04f\n", data[j]);
        if(data[j]>data[j-1])
            return TAIL;
    }
    return compute(data,k);

}
