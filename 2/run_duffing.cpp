#include <classical_oscillators.h>
#include <sstream>


int main(){
    double a,b,d,g,w;
    double x0, dx0;
    int steps;
    double dt;
    std::string file_output;

    

    std::ifstream file("input/duffing.txt");
    file >> a >> b >> d >> g >> w;
    file >> x0 >> dx0;
    file >> steps;
    file >> dt;
    file >> file_output;  
    file.close();

    std::string path = "output/duffing/" + file_output + ".txt";

    duffing_oscillator dos(d,a,b,g,w);
    dos.run_simulation(steps,dt,x0,dx0,path);
}