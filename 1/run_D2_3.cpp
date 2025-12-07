#include <random_placer.h>

int main(){
    int n;
    double density;
    int dim;
    
    n=40;
    dim=2;
    density=0.5200;
    std::string filepath = "output/statistics/" + std::to_string(density) + ".txt";
    std::ofstream outfile(filepath);
    for (int i=1;i<1000;i++){
        random_box some_elements(n,density,dim);
    
        outfile << some_elements.file_output3() << std::endl;
    
    }
}