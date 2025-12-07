#include <integrations.h>

double trapez_method(int N, double a, double b, std::function<double(double)> func){
    double h=(b-a)/N;
    double integral= 0.5 * (func(a) + func(b));

    for (int i = 1; i < N; ++i) {
        double x = a + i * h;
        integral += func(x);
    }
    return integral*h;
}
double trapez_method(int N[2], double a[2], double b[2], std::function<double(double, double)> func){//with a,b arrays
    double ha=(b[0]-a[0])/N[0];
    double hb=(b[1]-a[1])/N[1];
    double integral=0;
    for (int i = 0; i <= N[0]; ++i) {
        double x = a[0] + i * ha;
        for (int j = 0; j <= N[1]; ++j) {
            double y = a[1] + j * hb;
            
            double weight = 1.0;
            if (i == 0 || i == N[0]) weight *= 0.5; 
            if (j == 0 || j == N[1]) weight *= 0.5; 

            integral += weight * func(x, y);
        }
    }
    return integral*ha*hb;
} 

double simpson_method(int N, double a, double b, std::function<double(double)> func){
    
}
double simpson_method(int N[2], double a[2], double b[2], std::function<double(double, double)> func);
double monte_carlo(int N, double a, double b, std::function<double(double)> func);
double monte_carlo(int N[2], double a[2], double b[2], std::function<double(double, double)> func);
