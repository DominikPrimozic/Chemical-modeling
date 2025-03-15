#include <random_placer.h>

int main(){
    int n;
    double density;
    int dim;
    
    n=40;
    dim=2;
    density=0.1;
    for (int i=1;i<250;i++){
        density=i*0.005;
        std::string filepath = "output/crystals/" + std::to_string(density) + ".txt";
        random_box some_elements(n,density,dim);
        some_elements.file_output2(filepath);
    }
    
    
}