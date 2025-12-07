#include <2D_reactor.h>

int main(){
    int na,nb,np;
    double densitya,densityb;

    na=20;
    nb=20;
    np=5;
    densitya=0.25;
    densityb=0.25;

    porous_reactor leggo(na,nb,np,densitya,densityb,1);
    leggo.run(100,0.1);
    
}