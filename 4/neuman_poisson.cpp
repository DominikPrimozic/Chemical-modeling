#include <heat_conduction.h>
#include <fstream>
#include <string>

int main(){

    std::ifstream inputFile("input/poisson_neuman/input.txt");
    
    if (!inputFile) {
        std::cerr << "Error opening file.\n";
        return 1;
    }

    double left, right, bottom, top,source;
    int N;
    std::string line;
    std::getline(inputFile, line);
    left = std::stod(line);
    std::getline(inputFile, line);
    right = std::stod(line);
    std::getline(inputFile, line);
    bottom = std::stod(line);
    std::getline(inputFile, line);
    top = std::stod(line);
    std::getline(inputFile, line);
    source = std::stod(line);
    std::getline(inputFile, line);
    N = std::stoi(line);
    inputFile.close();

    std::string path="output/poisson_neuman/" + std::to_string(bottom) +"_"+ std::to_string(left)  +".txt";

    rectangular_rod rr(N,left,right,bottom,top);
    rr.solver_poisson_neuman(1.0/N,source);
    rr.print_to_file(path);


}

/////fuck meee, wron order for matrix