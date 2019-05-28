#include "filter.h"
#include "global.h"

void filter(double * data, double * filtered_data,int n){
    flip(data,filtered_data,n);

    for(int i; i<FILTER_AMOUNT;i++){
        average(filtered_data,filtered_data,n);
    }

}
void flip(double * data, double * flipped, int n){
    for(int i = 0;i<n;i++){
        flipped[i] = data[n-i-1];
    }
    return flipped;
}
void average(double * data, double *av, int n){
    for(int i=0;i<n;i++){
        int j=i;
        int k=0;
        double value=0;
        while(j>=0 && k<WINDOW){
            value+=data[j];
            j--;
            k++;
        }
        av[i]=value/k;
        printf("i=%d, orig=%lf, value=%lf, av=%lf,window=%d\n",i,data[i],value,av[i],(k));
    }
}

int test_filter(){
    test_flip();
    test_average();
    double v1[5] = {1, 2, 3, 4, 5};
}

int test_flip(){
    double vettore1[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double vettore2[8] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    printf("---  Test flip  ---\n");
    printf("\n---  v1  ---\n");
    for(int i=0;i<10;i++) {
        printf("%.0lf ", vettore1[i]);
    }
    double vettore1f[10];
    flip(vettore1,vettore1f,10);

    for(int i=0;i<10;i++) {
        printf("%.0lf ", vettore1f[i]);
    }
    printf("\n---  v2  ---\n");
    for(int i=0;i<9;i++) {
        printf("%.0lf ", vettore2[i]);
    }
    printf("\n");
    double vettore2f[9];
    flip(vettore2,vettore2f,9);
    for(int i=0;i<9;i++) {
        printf("%.0lf ", vettore2f[i]);
    }
    printf("\n\n");
    return 0;
}
int test_average(){
    double vettore[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    printf("---  Test average  ---\n");
    for(int i=0;i<10;i++) {
        printf("%lf ", vettore[i]);
    }
    printf("\n");
    double vettoreAV[10];
    average(vettore, vettoreAV,10);
    for(int i=0;i<10;i++) {
        printf("%lf ", vettoreAV[i]);
    }
    printf("\n\n");
}
