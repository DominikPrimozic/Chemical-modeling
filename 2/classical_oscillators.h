#ifndef oscillators
#define oscillators

#include <cmath>
#include <string>
#include <fstream>


class duffing_oscillator{
    private:
        double d,a,b,g,w;

        double du(double u, double v, double t);
        double dv(double u, double v, double t);
    public:
        void run_simulation(int steps, double dt, double x0, double dx0,  std::string path);
        duffing_oscillator(double d, double a, double b, double g, double w);

};

#endif // oscillators

