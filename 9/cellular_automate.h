#ifndef cellular_automata
#define cellular_automata


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <random>


inline int idx(int x, int y, int Nx) {
    return x + Nx * y;
}

struct grid_vector {
    int Nx, Ny;
    std::vector<int> g;
    grid_vector(int Nx, int Ny);
    int& operator()(int x, int y);
    int operator()(int x, int y) const;
};

class BZ_reaction{
    private:
        grid_vector current, next;

    public:
        BZ_reaction(int Nx, int Ny);
        void set_state(int x, int y, int value);
        void random_initialize(double p_excited = 0.01, double p_refractory = 0.05, int R = 10);
        std::vector<int> count_excited_neighbors(int x, int y, int R);
        void ca_step(int R, int k1, int k2, int g);
        void write_frame(int frame_number, int R);
        void run(int steps, int R, int k1, int k2, int g);

};


#endif // lattice_boltzman_method