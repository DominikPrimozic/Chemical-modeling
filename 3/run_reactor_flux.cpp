#include <diffusion_reaction.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


std::string removeComments(const std::string& line) {
    size_t pos = line.find("//");
    return (pos != std::string::npos) ? line.substr(0, pos) : line;
}


int main(){

    double k; //reaction constant
    double c0;
    double D; //diffusion coefficient
    double H; //spatial length
    double rxn_order; //reaction order
    double l; //thermal conductivity
    double T0,TH;
    double q; //first derivative of T
    double R;
    double Ea;
    double Ed;
    double heat; 
    double dx;
    double h;
    double s0;

    std::ifstream file("input/input.txt");
    std::vector<double*> variables = { &k, &c0, &D, &H, &rxn_order, &l, &T0, &TH, &q, 
        &Ea, &Ed, &heat, &dx, &h, &s0 };

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

    temperature_diffusion_reactor drq(D,k,rxn_order,H,c0,l,heat,T0,q,Ea,Ed);
    drq.solver_flux(dx,s0,h);

}