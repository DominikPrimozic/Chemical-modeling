#include <heat_conduction.h>
#include <fstream>
#include <string>

int main(){

    std::ifstream inputFile("input/neuman/input.txt");
    
    if (!inputFile) {
        std::cerr << "Error opening file.\n";
        return 1;
    }

    double left, right, bottom, top;
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
    N = std::stoi(line);
    inputFile.close();

    std::string path="output/neuman/" + std::to_string(bottom) +"_"+ std::to_string(left)  +".txt";

    rectangular_rod rr(N,left,right,bottom,top);
    rr.solver_laplace_neuman(1);
    rr.print_to_file(path);


}

/////fuck meee, wron order for matrix