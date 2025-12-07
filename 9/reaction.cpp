#include <cellular_automate.h>

int main() {
    const int Nx = 1000;
    const int Ny = 1000;
    const int R = 50;       
    const int steps = 5000;  
    BZ_reaction bz(Nx, Ny);

    
    bz.random_initialize(0.05, 0.2, R);

    bz.run(steps, R,2,3,25);

    std::cout << "Simulation complete. PGM files written.\n";
    return 0;
}