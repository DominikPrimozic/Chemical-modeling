#ifndef diffusivereaction_MM_kinetics
#define diffusivereaction_MM_kinetics

#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

class MMdr{
    private:
        double c0; //starting
        double c_0,c_h; //bounadary
        double D;
        double v_max, K_m;
        double h; //height
        double time;

        //dimensionless
        double Da,Pe;
        std::vector<double> z;


        double MM_reaction(double z_i);

    public:
        MMdr(double c0, double c_0, double c_h, double D, double v_max, double K_m, double h, double time);
        void solver( double dt, double dy);
        void print_solver(std::string path, double dt, double dy);
        void reset();

};


#endif // diffusivereaction_MM_kinetics

