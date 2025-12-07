#ifndef harmonic
#define harmonic

#include <complex>
#include <vector>
#include <cmath>
#include <eigensolver.h>
#include <fstream>
#include <iostream>
#include <string>


# define M_PI           3.14159265358979323846

class quantum_harmonic{
    private:
        matrix H0;
        matrix ap,am;

    public:
        quantum_harmonic(int num);
        matrix perturbation_3(double alpha);
        void solve_states(double alpha);

        void get_states(double alpha, double L, int N);

};

#endif // harmonic

