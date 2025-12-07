#ifndef MC_variational
#define MC_variational

#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <omp.h>


#define M_PI           3.14159265358979323846



class electron_in_potential{
    private:
        std::mt19937 gen;


        double trial_wavefunction(double r, double a);
        double local_energy(double r, double a); //for square potnetial
        bool metropolis(double& r, double a, double delta);


    public:
        double expected_r;
        double expected_r2;

        electron_in_potential();
        double run(double a, int steps, double delta, int discard=10000, int sample_interval=100, int tune=100);
        double optimize_alpha(double a_min, double a_max,double delta=1.0, double tol = 1e-3, int steps = 100000,int discard=10000, int sample=100, int tune=100); 
};


#endif // MC_variational

