#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
*/

/*
I use the Numerov method to solve the radial Schrodinger equation with the 3D Hamiltonian H = -1/2 d^2/dr^2 + 1/2 r^2
I look for energy levels in the interval [Emin, Emax] and I write on files the corresponding wave functions.
*/

double rmax = 8;
int n = 100000;
double h = rmax / n;
int counter = 0;
int l = 0;

/*
k_square corresponding to the 3D Hamiltonian H = -1/2 d^2/dr^2 + 1/2 r^2
*/
double k_square(double E, double r){
    return 2*E - pow(r,2) - l*(l+1)/pow(r,2);
}

/*
I calculate the wave function with the Numerov algorithm from 0 to rmax.
As initial condition I choose y(0) = 0 and y(h) = h^(l+1).
*/
double numerov(double E){
    double k_square0 = 0;
    double k_square1 = k_square(E, h);
    double k_square2 = k_square(E, 2*h);
    double y0 = 0;
    double y1 = pow(h, l+1);
    double y2 = y1 * (2-5*pow(h,2)*k_square1/6) / (1+pow(h,2)*k_square2/12);
    for(int i = 3; i <= n; i++){
        k_square0 = k_square1;
        k_square1 = k_square2;
        k_square2 = k_square(E, i*h);
        y0 = y1;
        y1 = y2;
        y2 = (y1 * (2-5*pow(h,2)*k_square1/6) - y0 * (1+pow(h,2)*k_square0/12))/ (1+pow(h,2)*k_square2/12);
    }
    return y2;
}

/*
The same as numerov but I normalize the wave function and save it to file.
Notice that if rmax is too big the function at rmax goes bananas.
Anyway the energy is more accurate with a big rmax.
*/
double numerovSave(double E, char * fileName){
    double k_square0 = 0;
    double k_square1 = k_square(E, h);
    double k_square2 = k_square(E, 2*h);
    double y[n+1];
    y[0] = 0;
    y[1] = pow(h,l+1);
    y[2] = y[1] * (2-5*pow(h,2)*k_square1/6) / (1+pow(h,2)*k_square2/12);
    for(int i = 2; i < n; i++){
        k_square0 = k_square1;
        k_square1 = k_square2;
        k_square2 = k_square(E, (i+1)*h);
        y[i+1] = (y[i] * (2-5*pow(h,2)*k_square1/6) - y[i-1] * (1+pow(h,2)*k_square0/12))/ (1+pow(h,2)*k_square2/12);
    }

    // I calculate the norm with the trapeziudal rule
    double norm = (pow(y[0],2) + pow(y[n],2)) / 2;
    for(int i = 1; i < n; i++){
        norm += pow(y[i],2);
    }
    norm = sqrt(norm * h);

    // I write on file
    FILE * file = fopen(fileName, "w");
    for(int i = 0; i <= n; i++){
        fprintf(file, "%.9f;%.9f\n", i*h, y[i]/norm);
    }
    fclose(file);

    return y[n];
}

/*
I use the secant method to find the zero of a function that changes sign in the interval [xmin, xmax].
The algorithm stops when I reach a certain precision.
In the end I print the obtained result.
*/
void secantMethod(double xmin, double xmax, double f_xmin, double f_xmax){
    double x;
    double precision = 0.000000001;
    while(fabs((xmax - xmin)/(xmax + xmin)) >= precision){
        x = xmax;
        xmax -= (xmax - xmin) / (f_xmax - f_xmin) * f_xmax;
        xmin = x;
        f_xmin = f_xmax;
        f_xmax = numerov(xmax);
    }
    char str [10];
    sprintf(str, "E%d_%d.txt\0", counter, l);
    numerovSave((xmax+xmin)/2, str);
    printf("E(n = %d, l = %d) = %.15f +- %.15f\n", counter, l, (xmax+xmin)/2, fabs(xmax-xmin)/2);
    counter ++;
}

/*
I divide the interval [xmin, xmax] in subintervals equally spaced dx.
For every subinterval where the function changes sign I apply the secantMethod.
*/
void zeroFinder(double xmin, double xmax, double dx){
    double f_xmin = numerov(xmin);
    double f_xmax;
    for(int i = 1; i <= (int)((xmax-xmin)/dx); i++){
        f_xmax = numerov(xmin + i*dx);
        if(f_xmin * f_xmax < 0){
            secantMethod(xmin+(i-1)*dx, xmin+i*dx, f_xmin, f_xmax);
        }
        f_xmin = f_xmax;
    }
}


int main(){
    double Emin = 1/(4*atan(1));
    double Emax = 8;
    double deltaE = 0.1;
    for(l = 0; l < 3; l++){
        counter = 0;
        zeroFinder(Emin, Emax, deltaE);
    }

    return 0;
}
