#include <2D_reactor.h>

int main(){
    int na,nb;
    double densitya,densityb;

    na=40;
    nb=40;
    densitya=0.25;
    densityb=0.25;
    //cutoff radius changing
    

    
    for (int i=40;i<80;++i){
        std::string filepath = "output/concentration_changes/a_025/" + std::to_string(i) + ".txt";
        reactor leggo(na,i,densitya,densityb,1);
        leggo.run2(250,0.2, filepath);
    }
    
}