#include <random_placer.h>

int main(){
    int n;
    double density;
    int dim;

    std::cout<<"Number of atoms: ";
    std::cin>> n;
    std::cout<<"Number of dimensions: ";
    std::cin>> dim;
    std::cout<<"Density: ";
    std::cin>> density;
    std::string filepath = "output/placements/"+ std::to_string(dim)+ "." + std::to_string(n)+ "." +std::to_string(density) + ".txt";
    random_box some_elements(n,density,dim);
    some_elements.file_output2(filepath);
        
}    