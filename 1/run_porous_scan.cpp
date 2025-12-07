#include <2D_reactor.h>

int main(){
    int na,nb,np;
    double densitya,densityb;

    na=20;
    nb=20;
    np=5;
    densitya=0.2;
    densityb=0.2;
    

    
    for (int i=0;i<21;++i){
        std::string filepath = "output/porous/concentrations/" + std::to_string(i) + ".txt";
        porous_reactor leggo(na-i,nb-i,2*i,densitya,densityb,1); //i want to keep constant concnetration of a and b
        leggo.run2(250,0.2, filepath);
    }
    
}