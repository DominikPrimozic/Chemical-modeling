#ifndef MD_simulation
#define MD_simulation

#include <vector>
#include <random>
#include <omp.h>
#include <VectorOps.h>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>
namespace fs = std::filesystem;

class LJ_cluster{
    private:
        std::vector<std::vector<double>> r,v,f, f_temp;
        int N,dim;

        //system variables
        double T;
        std::vector<double> force_temp;
        std::string position, velocity;
        std::string pat;
        bool track_positions,track_velocities,log;

        void rescale_thermostat();
        double kinetic_energy();
        double measure_temperature(double KE);
        double potential_energy(double epsilon=1, double sigma=1);


    public:
        LJ_cluster( std::string pat, bool log, int dim=2);
        double distance2(const std::vector<double>& p1, const std::vector<double>& p2);
        double pair_force(double r2, double sigma=1, double epsilon=1);
        void update_force();
        void velocity_verlet(double dt);
        void write_positions_binary(const std::string& filename);
        void write_velocities_binary(const std::string& filename);
        void run(double setT, int equilibration_steps, int steps, double dt, int thermostat_interval, int sampling_interval);
        void write_metadata(const std::string& filename, int equilibration_steps,int steps, double dt);
        double total_energy();


        void steepest_descent_minimize(double dt, int max_steps);
        void random_initialize(double box_size);
};


#endif // MD_simulation

