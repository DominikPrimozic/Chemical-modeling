#include <Michaelis_Menten_reactor.h>
#include <iostream>
#include <fstream>
#include <sstream>

std::string removeComments(const std::string& line) {
    size_t pos = line.find("//");
    return (pos != std::string::npos) ? line.substr(0, pos) : line;
}

int main(){

    double c0, c_0, c_h,  D, v_max,  K_m,  h,  time;
    double Ny,Nt;

    std::ifstream file("input/input.txt");
    std::vector<double*> variables = { &c0, &c_0, &c_h, &D, &v_max, &K_m, &h, &time, &Ny, 
        &Nt};

    std::string line;
    size_t index = 0;

    // Read and parse the file
    while (std::getline(file, line) && index < variables.size()) {
    line = removeComments(line);  // Remove comments
    std::istringstream iss(line);
    if (iss >> *variables[index]) {
    index++;
    }
    }

    file.close();

    double dy=1.0/(static_cast<int>(Ny));
    double dt=1.0/(static_cast<int>(Nt));

    std::cout << "Nt: " << dy << std::endl;
    std::cout << "Nt: " << dt << std::endl;

    MMdr reactor(c0,c_0,c_h,D,v_max,K_m,h,time);

    std::string path="output/reactor/" + std::to_string(D) + "_" + std::to_string(v_max) + "_" + std::to_string(K_m) + "_"+ std::to_string(dy) + "_" + std::to_string(dt) + ".txt";
    reactor.print_solver(path, dt,dy);


}