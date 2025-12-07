#include <cluster_dynamics.h>
#include <iostream>
#include <fstream>
#include <string>

struct SimulationConfig {
    int dimensions;
    std::string outputPath;
    int keepLog;
    double T_target;
    int equilibrationSteps;
    int productionSteps;
    double timeStep;
    int thermostatFrequency;
    int sampleFrequency;

};

#include <sstream> // For std::istringstream

bool parseConfigFile(const std::string& filename, SimulationConfig& config) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open config file: " << filename << std::endl;
        return false;
    }

    auto get_next_value = [&file](auto& var) -> bool {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            iss >> var;
            if (!iss.fail()) {
                return true;
            }
        }
        return false;
    };

    bool success = true;
    success &= get_next_value(config.dimensions);
    success &= get_next_value(config.outputPath);
    success &= get_next_value(config.keepLog);
    success &= get_next_value(config.T_target);
    success &= get_next_value(config.equilibrationSteps);
    success &= get_next_value(config.productionSteps);
    success &= get_next_value(config.timeStep);
    success &= get_next_value(config.thermostatFrequency);
    success &= get_next_value(config.sampleFrequency);

    if (!success) {
        std::cerr << "Error while parsing config file." << std::endl;
        return false;
    }

    return true;
}


int main(){
    SimulationConfig config;
    if (parseConfigFile("input/cluster/input.txt", config)) {
        std::cout << "Configuration loaded successfully!" << std::endl;
        // ... you can print the rest if you want
    } else {
        std::cerr << "Failed to load configuration." << std::endl;
    }


    LJ_cluster lj(config.outputPath,config.keepLog,config.dimensions);
    lj.run(config.T_target,    // T_target
        config.equilibrationSteps,    // equil steps
        config.productionSteps,   // production steps
        config.timeStep,  // dt
        config.thermostatFrequency,     // thermostat every 10 steps
        config.sampleFrequency);     // sample every 10 steps
    
}