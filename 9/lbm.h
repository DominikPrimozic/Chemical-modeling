#ifndef lattice_boltzman_method
#define lattice_boltzman_method

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdexcept>

inline int idx(int x, int y, int Nx) {
    return x + Nx * y;
}

struct lattice_velocities {
    int Nx, Ny;
    std::vector<double> f[9];
    lattice_velocities(int Nx, int Ny);
    double& operator()(int i, int x, int y);
    double operator()(int i, int x, int y) const;
    double& operator()(int i, int idx);
    double operator()(int i, int idx) const;
};

struct grid_vector {
    int Nx, Ny;
    std::vector<double> g;
    grid_vector(int Nx, int Ny);
    double& operator()(int x, int y);
    double operator()(int x, int y) const;
};

class D2Q9_2D {
private:
    lattice_velocities f, f0;
    std::vector<double> w;
    std::vector<int> cx, cy;
    std::vector<int> opp;
    double gravity, dt;
    grid_vector obstacle_mask;

public:
    int Nx, Ny;
    double tau; // relaxation time
    double nu;  // viscosity
    double dx;

    grid_vector rho, vx, vy;

    D2Q9_2D(int Nx, int Ny, double Re, double U, double L);
    void equilibrium_f();
    void collision();
    void stream();
    void bounce_boundaries();
    void macroscopic();
    void step();
    void run();
};

void write_field_to_file(const grid_vector& field, const std::string& filename, int step);
void print_grid_vector(const grid_vector& rho);


#endif // lattice_boltzman_method

 