#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
// #include <cmath>
#include <math.h>
// #include <cmath>

/* our custom libraries */
#include "save_data.h" // Include the custom header file
// #include "RngStream.h"


/* functions declaration for exercise A.3 */
double f_PDF(double x);
double g_PDF(double x, double p);
double g(double x);
double inverse_cubic(double a, double b, double u);


int main(){
    double a = 0.2;
    double b = 0.8;
    double c = 0.7;
    double d = 0.9;

    /*
    =========================================================
    = Sampling random points within D - dimensional domains =
    =========================================================
    */

    double Area = (b-a) * (d-c);
    double A_estimate = 0;
    double A_disk = 0;
    int R = 1;
    double Area_d = M_PI*(R*R);
    fprintf(stdout,"Rectangle exercise: \n\n");
    
    int N = 1000000; // Number of throws

    // Set up the rnd generator
    int iseed=0;
    srand(iseed);
    
    double x, y = 0;

    for (int i = 0; i < N; i++){
        x = drand48();
        y = drand48();
        if ((x < b & x > a) & (y < d & y > c)){
            A_estimate += 1;
        }
        if (sqrt(x*x + y*y) <= R){
            A_disk += 1;
        }
    }
    fprintf(stdout,"Real area = %e\n", Area);
    fprintf(stdout, "The estimated area is : %e\n\n", A_estimate/N);
    fprintf(stdout,"Real disk area = %e\n", Area_d);
    fprintf(stdout, "The estimated area is : %e\n\n", 4*A_disk/N);
    
    /*
    ====================
    = Inversion method =
    ====================
    */

    fprintf(stdout, "\n- Exercise 1\n");

    int n = 3; // case n = 3, 4
    // from normalization 
    c = n+1;
    double *vec = (double*)malloc(N * sizeof(double));
    double u = 0;

    for (int i = 0; i < N; i++){
        u = drand48(); // sample from the cumulative
        // The inverse of the cumulative is x = u^(1/(n-1))
        vec[i] = pow(u, 1.0/(n+1));
    }
    // save results in a file
    char s[100] = "output_ex_2_2_1.txt";
    saveResultsToFile(vec, N, s);  
    // Results in "ex_n_2_2_1.png"

    // free memory
    // free(vec);
    fprintf(stdout, "Done! Output file saved!\n");

    fprintf(stdout, "\n- Exercise 2:\n");
    // Now, I want to estimate using MonteCarlo the normalization constant c
    double theor_c = 3./8.;
    double integral = 0;
    double M = 4;
    double m = 0;
    a = 0; 
    b = 2;
    for (int i = 0; i < N; i++){
        // Hit and miss method:
        u = a+(b-a)*drand48(); // x
        y = (M-m)*drand48(); // (M-m)*u
        if (y+m < u*u){
            integral += 1;
        }
    }
    integral = (integral / N) * (b-a)*(M-m) + m*(b-a);
    // Since the distribution must be normalized:
    c = 1/integral;
    fprintf(stdout, "Theoretical c = %e\nEstimated c = %e\n", theor_c, c);
    // Use theoretical c:
    c = theor_c;
    // Sample from the distribution cx^2
    for (int i = 0; i < N; i++){
        u = drand48();
        vec[i] = cbrt(3.*u/c);
    }
    // save results in a file
    char st[100] = "output_ex_2_2_2.txt";
    saveResultsToFile(vec, N, st); 
    // Results in "ex_n_2_2_2.png"

    // free memory
    // free(vec);
    fprintf(stdout, "Done! Output file saved!\n");

    // Additional exercises
    // 1.
    double mu = 1.0;
    for (int i = 0; i < N; i++){
        u = drand48(); // sample from the cumulative
        // The inverse of the cumulative is x = -log(1-u)/mu
        vec[i] = -log(1-u)/mu;
    }
    char add_1[100] = "additional_1.txt";
    saveResultsToFile(vec, N, add_1); 
    // 2.
    for (int i = 0; i < N; i++){
        u = drand48(); // sample from the cumulative
        // The inverse of the cumulative is x = sqrt(-log(1-u))
        vec[i] = sqrt(-log(1-u));
    }
    char add_2[100] = "additional_2.txt";
    saveResultsToFile(vec, N, add_2); 

    // 3.
    a = 1.;
    n = 3;
    b = - pow(a, 1.-n)/(1.-n);
    for (int i = 0; i < N; i++){
        u = drand48(); // sample from the cumulative
        // The inverse of the cumulative is x = {[(u + a^{1-n}/(b(1-n)))* b(1-n)]^{1/(1-n)} - a}/b
        x = (pow(u*b*(1-n)+pow(a, 1-n), 1./(1-n))-a)/b;
        // if (x <= 4){
        vec[i] = x;
        // }
    }
    char add_3[100] = "additional_3.txt";
    saveResultsToFile(vec, N, add_3); 

    /*
    ==============================================
    = Sampling via transformation of coordinates =
    ==============================================
    */

    fprintf(stdout, "\n- Exercise 1\n");
    fprintf(stdout, "Wrong disk sample\n");

    // Wrong sample

    M = 2;
    double **matrix = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        matrix[i] = (double *)malloc(M * sizeof(double));
    }

    for (int i = 0; i < N; i++){
        matrix[i][0] = drand48();
        matrix[i][1] = 2*M_PI*drand48();
    }
    char s2[100] = "wrong.txt";
    saveMatrixToFile(matrix, N, M, s2);
    // In the plot 'wrong_uniform.png' we can see how the data plotted are not uniform on a circle
    // but there is an higher concentration of points around the center.

    fprintf(stdout, "Correct disk sample\n");

    /*
    The conceptual mistake of this algortihm is to not have considered the dipendence of r and θ
    In cartesian coordinates we can use this method to generate x = r*cos(θ) and y = r*sin(θ), 
    but θ = tan^-1 (y/x) depends on both x and y as r = sqrt(x^2 + y^2) and it is not fully independedent.
    
    In order to fix this, we correctly apply the change of variables:
    p(r, θ) = p(x, y)*r = r/(πR^2)
    and p(r) = integral{p(r, θ) dθ} = 2*r/(R^2), where r ∈ [0, R]
    Finally, p(θ|r) = p(r, θ)/p(r) = 1/(2π)
    
    In order to sample from p(r), the cumulative P(r) is given by: P(r) = (r/R)^2 = u 
                                                                    <==> r = sqrt(u * R^2) 
    Instead, the cumulative of p(θ|r) is: P(θ|r) = θ/(2π) = u' 
                                                                    <==>  θ = 2πu'
    where u, u' ∈ [0, 1]                                                             
    */

    for (int i = 0; i < N; i++){
        matrix[i][0] = sqrt(drand48() * R*R); // r sampled from p(r)
        matrix[i][1] = 2*M_PI*drand48(); // θ sampled from p(θ|r)
    }
    char s3[100] = "correct_disk_sampling.txt";
    saveMatrixToFile(matrix, N, M, s3);
    // Results in "correct_disk_sampling.png"

    fprintf(stdout, "\n- Exercise 2: Gaussian sample\n");

    /*
    p(r, θ) dr dθ = 1/(2π) exp(-(r^2)/2) r dr dθ = p(r) p(θ) dr dθ
    Therefore, we can sample r ~ p(r) and θ ~ p(θ)

    In order to draw r ~ p(r), the cumulative P(r) is given by:
    P(r) =  integral(r' exp(-(r'^2)/2) dr') from r' = 0 to r' = r
         = 1 - exp(-(r^2)/2) = u 
         <==> r = sqrt(-2*log(1-u))
    To draw θ ~ p(θ), the cumulative P(θ) is of course:
    P(θ) = integral(dθ'/(2π)) from θ' = 0 to θ' = θ
         = θ/(2π) = u' 
         <==> θ = 2πu'
    with u, u' ∈ [0, 1]    
    */

    for (int i = 0; i < N; i++){
        matrix[i][0] = sqrt(-2*log(1-drand48())); // r sampled from p(r)
        matrix[i][1] = 2*M_PI*drand48(); // θ sampled from p(θ)
    }
    char s4[100] = "gaussian_sample.txt";
    saveMatrixToFile(matrix, N, M, s4);
    // Results in "gaussian_sampling.png"

    // free memory
    free(matrix);

    /*
    The above can be easily extended to the most general case of N(μ, σ^2):
    starting from p(r, θ) above defined, the new samples x' = σx + μ = σ rcos(θ) + μ  
    and y' = σ rsin(θ) + μ, where r, θ ~ p(r, θ). Therefore, having sampled from N(0, 1) we can get any Gaussian 
    of mean μ and variance σ^2.
    */

    fprintf(stdout, "\n- Exercise 2: Rejection method\n");
    // In order to be normalized, the g(x) function has constant value A = 2p/(2p^2 + 1)
    fprintf(stdout, "\n- Evaluate best p\n");
    int j = 0;
    double p = 0.;
    double A = 0.;
    // create a matrix
    double **mat = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        mat[i] = (double *)malloc(M * sizeof(double));
    }
    double xi = 0.; // for the Von Neumann method 
    double add = 0.;
    int k = 0;
    char s5[100] = "rejection_method.txt";
    double *eff = (double*)malloc(30 * sizeof(double));
    while (j < 30){
        add = 0;
        p += 0.1;
        A = (2.*p)/(2*p*p + 1);
        // The c constant can have multiple possible values 
        // For x ∈ [0, p], to achieve f(x) <= c*g(x) I set c >= sqrt(2/π)/A 
        c = sqrt(2/M_PI)/A;
        // The cumulative of g(x) if the function G(x) = Ax for 0 <= x <= p and G(x) = 1- A/(2p) * exp(p^2 - x^2) for x > p 
        // Applying the inverse method to generate a sample X from g(x): x = u/A if 0 <= u <= pA and 
        // x = sqrt(p^2 - log(2(1-u)*p/A)) for u > pA

        if (p == 0.2){
            while (k < N){
                u = drand48();
                xi = drand48();
                if (u <= p*A){
                    x = u/A;
                }
                else{
                    x = sqrt(p*p - log(2*p*(1-u)/A));
                }
                if (xi < f_PDF(x)/(c*g_PDF(x, p))){
                    mat[k][0] = x;
                    mat[k][1] = xi;
                    k += 1;
                }
            }
            saveMatrixToFile(mat, N, M, s5);
            // Results in "rejection_method.png" 
        }
        for (int k = 0; k < N; k++){
            u = drand48();
            xi = drand48();
            if (u <= p*A){
                x = u/A;
            }
            else{
                x = sqrt(p*p - log(2*p*(1-u)/A));
            }
            if (xi < f_PDF(x)/(c*g_PDF(x, p))){
                add += 1;
            }
        }
        eff[j] = add;
        j += 1;
    }
    char s_p[100] = "find_best_p.txt";
    saveResultsToFile(eff, 30, s_p);
    
    free(eff);

    /*
    =======================
    = Importance sampling =
    =======================
    */

    fprintf(stdout, "\n- Exercise 1 - Crude approach\n");
    b = 10;
    integral = 0;
    a = 0;
    for (int i = 0; i < N; i++){
        x = (b-a)*drand48() + a;
        integral += exp(-x*x)*g(x);  
    }
    integral = integral*(b-a)/N;
    fprintf(stdout, "Estimated integral via crude MC = %e with g(x) = x^2/(x^2+1)\n", integral); // Theoretical value 0.214580214629...
    
    fprintf(stdout, "\n- Exercise 1 - importance sampling\n");
   
    /*
    W(x) is the f(x) from which we have sampled in the previous exercise: since now the requirement is that W(x)
    is a PDF, it must me normalized in the [0, +∞] interval. So W(x) -> A*W(x) where A = sqrt(2)
    This means that after evaluating g(x_i) and summing for i going from 0 to N, we have to multiply not only 
    for sqrt(π/2) but also for 1/A, so for sqrt(π)/2 = sqrt(π/4) 
    */

    integral = 0;
    for (int i = 0; i < N; i++){
        x = mat[i][0];
        integral += g(x);  
    }
    integral = integral*sqrt(M_PI/4)/N;
    fprintf(stdout, "Estimated integral via importance sampling = %e with g(x) = x^2/(x^2+1)\n", integral);


    fprintf(stdout, "\n- Exercise 2 - cos(x) integral\n");
    double theor_integral = 1;
    
    /*
    In this very last exercise the goal is to apply the same method as before, but with f(x) = cos(x) and
    g(x) = a + bx^2, which has cumulative G(x) = ax + bx^3/3 for x ∈ [0, π/2]. First thing to do is to invert 
    a cubic. The way we choose a, b is given by the fact that we want h(x) = f(x)/g(x) ~ const. and the constraint
    that g(x) in [0, π/2] must be normalized. 
    Therefore, a and b are linked and the must approximate the f(x) = cos(x) for x ∈ [0, π/2] .
    This is achieved by requiring that g(0) = 1 => a = 1 and g(π/2) = 0 => b = -4/π^2.
    Of course, given b, a ~ 1 since the normalization constraint must be satisfied.
    Then, the sum over h(x_i) (i ∈ [1, k]) is performed and the final integral is given by dividing by k:
    k represents the number of iterations 1 <= k <= N in order to have an accuracy ~ 1%.
    */
 
    b = -4/(M_PI*M_PI);
    fprintf(stdout, "b = %e\n", b);
    a = (1-b*pow(M_PI, 3)/24)*2/M_PI; // g(x) must be normalized
    fprintf(stdout, "a = %e\n", a);
    
    integral = 0;
    k = 0;
    double accuracy = 0;
    while (k < N){
        u = drand48();
        x = inverse_cubic(a, b, u); // x is sampled according to g(x)
        if (x > 0){
            integral += cos(x)/(a+b*x*x);
            k += 1;
        }
        accuracy = fabs(integral/k-theor_integral)/(theor_integral);
        if (accuracy <= 0.01){
            break;
        }
    }
    integral /= k;
    fprintf(stdout, "The integral is theoretically = %e\n", theor_integral);
    fprintf(stdout, "The estimate via the importance sampling method is = %e \n", integral);
    fprintf(stdout, "With an accuracy of = %e, after %d iterations! \n\n", accuracy, k); // ~ 3 iterations to compute the integral


    // free memory
    free(vec);
    free(mat);
    return 0;


}

/* functions definition for exercise 2.3*/
double f_PDF(double x){
    return sqrt(2/M_PI) * exp(-x*x);
}
double g_PDF(double x, double p){
    double A = (2.*p)/(2*(p*p) + 1);
    if ((x <= 0) & (x <= p)){
        return A;
    }
    return A/p * x * exp(p*p - x*x);
}

/* slowly varying function: Exercise 1 importance sampling */
double g(double x){
    return x*x/(x*x+1);
    // return log(x+1);
    // return 1;
}

/* inverse cubic function */
double inverse_cubic(double a, double b, double u){
    // convert to a depressed cubic equation:
    double p = 3*a/b;
    double q = -3*u/b;
    double D = (q*q)/4 + (p*p*p)/27;
    double out = 0;
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
        out = out2; // hopefully this is the right one
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