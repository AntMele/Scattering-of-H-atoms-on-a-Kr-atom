#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double rmax = 15;
int n = 250000;
double h = rmax / n;
int l = 0;
double beta = 0.03518;
double eta = 0.3;
double k = sqrt(eta / beta); // E = hbar^2 * k^2 / 2m but actually in unit of sigma it is k*sigma because k*r=k*sigma*x


double sinc(double x){
    return sin(x)/x;
}

double cosc(double x){
    return cos(x)/x;
}

double jbesseln(double x, int n){
    double y;
    if(n < 1){
        if(n == 0){
            return sinc(x);
        }else{
            return cosc(x);
        }
    }else{
        return ((2*(n-1)+1)/x)*jbesseln(x,n-1) - jbesseln(x,n-2);
    }
}

double nbesseln(double x, int n){
    if(n < 1){
        if(n == 0){
            return -cosc(x);
        }else{
            return sinc(x);
        }
    }else{
        return ((2*(n-1)+1)/x)*nbesseln(x,n-1) - nbesseln(x,n-2);
    }
}


double k_square(double x){
    return (eta - 4*(1/pow(x,12) - 1/pow(x,6))) / beta - l*(l+1)/pow(x,2);
}

void numerov(double* y, double xlow){
    y[0] = 0;
    int start = (int)(xlow/h);
    for(int i = 1; i <= start+2 ; i++)
        y[i] = exp(-2/(5*sqrt(beta)*pow(i*h,5)));
    double k_square0 = k_square(start*h);
    double k_square1 = k_square((start+1)*h);
    double k_square2 = k_square((start+2)*h);
    for(int i = start+2; i <= n; i++){
        k_square0 = k_square1;
        k_square1 = k_square2;
        k_square2 = k_square((i+1)*h);
        y[i+1] = (y[i] * (2-5*pow(h,2)*k_square1/6) - y[i-1] * (1+pow(h,2)*k_square0/12)) / (1+pow(h,2)*k_square2/12);
    }

    /*FILE * file = fopen("l00.txt", "w");
    for(int i = 0; i <= n; i++){
        fprintf(file, "%.9f;%.9f\n", i*h, y[i]);
    }
    fclose(file);*/
}

double pshift(double* y, double x1, double x2){
    double key = x2 * y[(int)(x1/h)] / (x1 * y[(int)(x2/h)]);
    return atan((key*jbesseln(k*x2,l) - jbesseln(k*x1,l)) / (key*nbesseln(k*x2,l) - nbesseln(k*x1,l)));
}

/*
int main(){
    double x1 = 14;
    double x2 = 14.2;
    double xlow = 0.4;
    FILE * file = fopen("psn.txt", "w");
    for(int i = 0; i <= 100; i++){
        n = round(750 * pow(1.06, i));
        h = rmax / n;
        fprintf(file, "%d", n);
        for(l = 0; l < 7; l++){
            numerov(y, xlow);
            fprintf(file, ";%.15f", pshift(y, x1, x2));
        }
        fprintf(file, "\n");
        printf("%d\n", n);
    }
    return 0;
}*/

int main(){
    double y[n+1];
    double x1 = 14;
    double x2 = 14.2;
    double xlow = 0.4;

    for(l = 0; l < 15; l++){
        numerov(y, xlow);
        printf("%.15f\n", pshift(y, x1, x2));
    }

    return 0;
}
