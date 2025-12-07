#include <monte_carlo.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>



int main() {
    int num_particles;
    double density, epsilon, sigma, cutoff, temperature;
    int dimensions, equilibration_steps, sampling_steps, tuning_steps;
    int equil_steps_check, sampling_interval;
    std::string energies_path, config_path, log_path;
    bool save_configs,keep_log;

    simulation_parameters sp;

    std::ifstream file("input/LJ.txt");

    if (!file.is_open()) {
        std::cerr << "Failed to open file!" << std::endl;
        return 1;
    }

    std::string line;
    int line_number = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        switch (line_number) {
            case 0: iss >> num_particles; break;
            case 1: iss >> density; break;
            case 2: iss >> dimensions; break;
            case 3: iss >> temperature; break;
            case 4: iss >> epsilon; break;
            case 5: iss >> sigma; break;
            case 6: iss >> cutoff; break;
            case 7: iss >> equilibration_steps; break;
            case 8: iss >> sampling_steps; break;
            case 9: iss >> tuning_steps; break;
            case 10: iss >> equil_steps_check; break;
            case 11: iss >> sampling_interval; break;
            case 12: iss >> energies_path; break;
            case 13: iss >> config_path; break;
            case 14: {
                int tmp;
                iss >> tmp;
                save_configs = (tmp == 1);
                break;
            }
            case 15: iss >> log_path; break;
            case 16: {
                int tmp;
                iss >> tmp;
                keep_log = (tmp == 1);
                break;
            }
            default: break;
        }
        ++line_number;
    }

    file.close();

    // Assign to struct
    sp.equil = equilibration_steps;
    sp.steps = sampling_steps;
    sp.tune = tuning_steps;
    sp.sample_interval = sampling_interval;
    sp.equil_steps_check = equil_steps_check;
    sp.energies = energies_path;
    sp.configurations = config_path;
    sp.save_configs = save_configs;
    sp.log_path = log_path;
    sp.keep_log = keep_log;


    Lennard_Jones lj(num_particles,density,temperature,epsilon,sigma,dimensions,cutoff);
    lj.run_simulation(sp);
    return 0;
}