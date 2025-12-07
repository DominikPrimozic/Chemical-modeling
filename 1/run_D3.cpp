#include <2D_reactor.h>

int main(){
    int na,nb;
    double densitya,densityb;

    std::cout<<"Number of atoms A: ";
    std::cin>> na;
    std::cout<<"Number of atoms B: ";
    std::cin>> nb;
    std::cout<<"Density A: ";
    std::cin>> densitya;
    std::cout<<"Density B: ";
    std::cin>> densityb;

    reactor leggo(na,nb,densitya,densityb,1);
    leggo.run(1000,0.1);
    
}