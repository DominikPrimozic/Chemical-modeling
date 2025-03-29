#ifndef diffusivereaction
#define diffusivereaction

#include <cmath>
#include <string>
#include <fstream>
#include <linear_solvers.h>
//will make my own solver but eigen is probably better, will try to add eigen functionality into this

class diffusion_reactor{
    private:
        double k; //reaction constant
        double c0;
        double D; //diffusion coefficient
        double H; //spatial length
        double rxn_order; //reaction order
        std::string outpath;

        double du(double x, double u, double v);
        double dv(double x, double u, double v);
        double RK4_solver(double dx, double s);
        double root_solver(double s, double h,double dx);
        double root_function(double s,double dx);

    public:

        diffusion_reactor(double D, double k, double a, double H, double c0);

        void final_RK4_solver(double dx, double s);


        void solver(double dx, double s0, double h);



};

class temperature_diffusion_reactor{
    private:
        double k; //reaction constant
        double c0;
        double D; //diffusion coefficient
        double H; //spatial length
        double rxn_order; //reaction order
        double l; //thermal conductivity
        double T0,TH;
        double q; //first derivative of T
        std::string outpath;

        double R;
        double Ea;
        double Ed;
        double heat; //reaction heat


        double reaction_constant(double w);
        double diffusion_constant(double w);

        double du(double x, double u, double v, double w, double z);
        double dv(double x, double u, double v, double w, double z);
        double dw(double x, double u, double v, double w, double z);
        double dz(double x, double u, double v, double w, double z);
        double RK4_solver_flux(double dx, double s);
        double root_solver_flux(double s, double h,double dx);
        double root_function_flux(double s,double dx);

        double temperature_profile(double x);
        double du(double x, double u, double v);
        double dv(double x, double u, double v);
        double RK4_solver(double dx, double s);
        double root_solver(double s, double h,double dx);
        double root_function(double s,double dx);

    public:

        //really since we have a reaction it produces heat that is a fucntion of concentracion, and we know T(0) and dT(0)
        temperature_diffusion_reactor(double D, double k, double a, double H, double c0, double l, double rheat, double T0, double q, double Ea, double Ed);
        void final_RK4_solver_flux(double dx, double s);
        void solver_flux(double dx, double s0, double h);

        //simplified version where c depends on T but T is only a function of x, so we first solve for T(x) and then for concnetration
        temperature_diffusion_reactor(double D, double k, double a, double H, double c0, double l,double T0, double TH,double Ea, double Ed); 
        void final_RK4_solver(double dx, double s);
        void solver(double dx, double s0, double h);

};

#endif // diffusivereaction

