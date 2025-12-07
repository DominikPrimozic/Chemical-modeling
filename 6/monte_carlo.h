#ifndef MC_simulation
#define MC_simulation

#include <vector>
#include <random_placer.h>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <filesystem>
#include <omp.h>


#define M_PI           3.14159265358979323846


struct simulation_parameters{
    int steps;
    int equil;
    int tune;
    double step_size;
    int sample_interval;
    int equil_steps_check;
    std::string energies;
    std::string configurations;
    std::string log_path;
    bool save_configs = true;
    bool keep_log = true;
};

class Lennard_Jones{
    private:
        int N;
        double density;
        double T;
        int dimensions;

        double e,s;
        double s2;
        double r_cutoff,r_cutoff2;

        std::vector<std::vector<double>> particles; //particle positions
        double L; //box size
        std::vector<double> trial_position; //for allocating

        double E;
        std::mt19937 gen;

        double distance2(std::vector<double>& p1, std::vector<double>& p2);
        double LJ_pair(double r2);
        double LJ_potential(int i, std::vector<double>& i_coords);
        double compute_energy();
        bool move(double step_size);
        bool metropolis(double dE);

    public:
        Lennard_Jones(int N, double density, double T, double e, double s, int dim, double rc);
        void run_simulation(simulation_parameters& sp);
        void equilibrate(simulation_parameters& sp);
        void sampling(simulation_parameters& sp);

        void sample_step(int step,  std::ofstream& out_energy,const simulation_parameters& sp);
        void write_particles_binary(const std::string& filename);
};

#endif // MC_simulation

