#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
// #include <cmath>
#include <math.h>
// #include <cmath>

double inverse_cubic(double a, double b, double u);


int main(){
    double b = -4/(M_PI*M_PI);
    double a = (1-b*pow(M_PI, 3)/24)*2/M_PI;
    double x = 0;
    fprintf(stdout, "a = %e \n", a);
    fprintf(stdout, "b = %e \n", b);
    for (int i = 0; i < 15; i++){
        double u = drand48();
        x = inverse_cubic(a, b, u); // x is sampled according to g(x)
    
        fprintf(stdout, "%d. u = %e and x = %e \n", i, u, x);

    }


    return 0;
}

/* inverse cubic function */
double inverse_cubic(double a, double b, double u){
    // convert to a depressed cubic equation:
    double p = 3*a/b;
    double q = -3*u/b;
    double D = (q*q)/4 + (p*p*p)/27;
    double out = 0;
    fprintf(stdout, "D = %e \n", D);
    if (D > 0){
        // In this case, there is only one real root, given by "out" below
        double u = cbrt(-q/2 + sqrt(D));
        double v = cbrt(-q/2 - sqrt(D));
        out = u+v;
    }
    else if (D < 0){
        /*
        In this case there 3 real distinct solutions, given by out1,
        out2, out3 below. The one that interests us is that in the
        inerval [0,1]. It is seen ("empirically") that is always the
        second one in the list below [there is perhaps more to search here]
        */
        double out1 = 2*sqrt(-p/3)*cos((acos(3*q/(2*p)*sqrt(-3/p)))/3);
        double out2 = 2*sqrt(-p/3)*cos((acos(3*q/(2*p)*sqrt(-3/p))/3-2*M_PI/3));
        double out3 = 2*sqrt(-p/3)*cos((acos(3*q/(2*p)*sqrt(-3/p))/3-4*M_PI/3));

        /*
        ToDo: implement a check just to be sure out2 is the good root 
        (in case this "empirical" truth turns out to stop working) 
        */
        out = out1; // hopefully this is the right one
    }
    else{
        /*
        In theory we always go from D>0 to D<0 by passing to a D=0
        boundary, where we have two real roots (and where the formulas
        above change again slightly). In practice, however, due to round-off errors,
        it seems we never hit this boundary but always pass "through" it 
        This D=0 scenario could still be implemented if needed, though.
        */
        double out1 = 3*q/p;
        double out2 = -3*q/(2*p);   // = out3
        if (out1 > 0){
            out = out1;
        }
        else{
            out = out2;
        }
    }
    return out; // return x sampled
}