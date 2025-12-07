#ifndef steady_state_conduction
#define steady_state_conduction

#include <cmath>
#include <string>
#include <fstream>
#include <linear_solvers.h>
#include <vector>


class rectangular_rod{
    private:
        matrix temperatures; //its row1, row2, row3,...
        
        int Nx,Ny; 
        int N; //assuming equal, this is a single dimension

        double x0,xn,y0,yn;

        matrix starting_temeperatures;

    public:
        rectangular_rod(int N, double x0, double xn, double y0, double yn);

        void boundary_values(double h);
        void boundary_values_neuman(double h);
        void boundary_values_robin(double h);
        void print_to_file(std::string path); //should just use matrix class and add to it access to data for GC_solver
        void solver_laplace(double h); //assumin dx=dy
        void solver_laplace_neuman(double h);

        void solver_poisson(double h,double source);
        void solver_poisson_neuman(double h, double source);

};


#endif // steady_state_conduction
