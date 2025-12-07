#include <distributions.h>
#include <sstream>
#include <iostream>


int main(){
    int N;
    double T,m;
    double maxV, minV;
    std::string file_output;

    

    std::ifstream file("input/maxwell.txt");
    file >> N;
    file >> m >> T;
    file >> minV >> maxV;
    file >> file_output;  
    file.close();

    m/=6.02e26;

    std::string path = "output/maxwell/" + file_output;

    maxwell_distribution md(N,path,T,m,minV,maxV);
}