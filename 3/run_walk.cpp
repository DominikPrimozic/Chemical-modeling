#include <random_walk.h>
#include <omp.h>

int main(){

    //walker letswalk(100,0.3);
    //letswalk.walk(100,0.5);
    //letswalk.multiple_walks(1000,100,0.5);

    int N=100;
    //double density;
    //only care about density chaninig
    /*
    #pragma omp parallel for private(density)
    for (int i=1;i<100;i++){
        density=i*0.006;
        walker letswalk(N,density);
        letswalk.multiple_walks(300,500,0.5);
    }*/


    //for second part i need some changes because i need to make 2D potential surface for changing both density and step size
    #pragma omp parallel for
    for (int i=1;i<50;i++){
        double density=i*0.012;
        walker letswalk(N, density);
        letswalk.update_paths(0.12);
        for (int j=1;j<50;j++){
            double step=j*0.04;
            letswalk.multiple_walks2(100, 250, step);
        }
        
    }

}