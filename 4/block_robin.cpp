#include <heat_conduction.h>
#include <fstream>
#include <string>

int main(){

    std::ifstream inputFile("input/convection/input.txt");
    
    if (!inputFile) {
        std::cerr << "Error opening file.\n";
        return 1;
    }

    double left, right, bottom, top,k,h;
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
    h = std::stod(line);
    std::getline(inputFile, line);
    k = std::stod(line);
    std::getline(inputFile, line);
    N = std::stoi(line);
    inputFile.close();

    std::string path="output/convection/" + std::to_string(bottom) +"_"+ std::to_string(left)  +".txt";

    rectangular_rod rr(N,left,right,bottom,top);
    rr.solver_laplace_robin(0.1, h/k,1,h/k*top);
    rr.print_to_file(path);


}

/////fuck meee, wron order for matrix