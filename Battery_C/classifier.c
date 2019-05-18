
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
    //test_classificatore();
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

}

int is_increasing(double *data, int n){
    for(int i = 1;i<n;i++){
        if(data[i-1]>data[i])
            return i;
    }
    return 0;
}




int is_decreasing(double *data, int n){
    for(int i = 1;i<n;i++){
        if(data[i-1]<data[i])
            return i;
    }
    return 0;
}

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
                return 1;
            }
        }
    }
    return 0;
}

int classificatore(double *data, int n){
    printf("-- debug classificatore --\n");
    int h = 0;
    int k = 0;
    int esito1 = 1;
    double m_array[n];
    while(esito1==1 && k+4<n){
        if(data[k]>data[k+1])
            esito1=0;
        double m = (data[k+4]-data[k])/5;
        m_array[h]=m;
        h++;
        k++;
    }
    int esito = compute(m_array,h);
    int j = 1;
    int esito2 = 1;
    printf("h: %d, esito: %d\n", h,esito);
    while(esito2==1 && k+j<n){
        if(data[k+j-1]<data[k+j])
            esito2=0;
        j++;
    }
    printf("J = %d, n = %d\n",j,n);
    if((j+k)==(n-1) && j>1){
        if (esito==1){
            return 2;
        } else {
        return 1;
        }
    } else {
        return 3;
    }
}
