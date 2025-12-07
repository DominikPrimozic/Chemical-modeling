#include <fstream>
#include <string>
#include <sstream>
#include <iostream>


struct reactor{
    double k,P,s,c0,T0;

    double du(double u, double t){
        return -k*(u-T0) + P/c0 - s*(u*u*u*u - T0*T0*T0*T0);
    }
};



int main(){

    double k,P,s,c0,T0 ;
    double x0;
    int steps;
    double dt;
    std::string file_output;

    reactor r;

    std::ifstream file("input/reactor.txt");
    file >> r.k>>r.P>>r.s>>r.c0>>r.T0;
    file >> x0;
    file >> steps;
    file >> dt;
    file >> file_output;  
    file.close();

    std::string path="output/reactor/" + file_output + ".txt";

    double t=0;
    double u=x0;
    std::ofstream outFile(path);
    outFile << t << "\t" << u << std::endl;
    for (size_t i=0;i<steps;++i){
        double u1= dt*r.du(u,t);
        double u2= dt*r.du(u+0.5*u1,t+0.5*dt);
        double u3= dt*r.du(u+0.5*u2,t+0.5*dt);
        double u4= dt*r.du(u+u3,t+dt);

        u+=1/6.0 * (u1+2*u2+2*u3+u4);
        t+=dt;

        outFile << t << "\t" << u << std::endl;
    }

    outFile.close();

}