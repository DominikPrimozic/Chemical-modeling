#ifndef integrals
#define integrals

#include <iostream>
#include <functional>
//lets do some overloading shananigans
double trapez_method(int N, double a, double b, std::function<double(double)> func);
double trapez_method(int N[2], double a[2], double b[2], std::function<double(double, double)> func); //with a,b arrays
double simpson_method(int N, double a, double b, std::function<double(double)> func);
double simpson_method(int N[2], double a[2], double b[2], std::function<double(double, double)> func);
double monte_carlo(int N, double a, double b, std::function<double(double)> func);
double monte_carlo(int N[2], double a[2], double b[2], std::function<double(double, double)> func);

#endif // integrals