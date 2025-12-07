#include <classical_oscillators.h>


double duffing_oscillator::du(double u, double v, double t){
    return v;
}

double duffing_oscillator::dv(double u, double v, double t){
    return g*std::cos(w*t) - d*v - a*u - b*u*u*u;
}

void duffing_oscillator::run_simulation(int steps, double dt, double x0, double dx0, std::string path ){
    double t=0;
    double u=x0, v=dx0;

    std::ofstream outFile(path);
    outFile << t << "\t" << u << "\t" << v << std::endl;
    for (size_t i=0;i<steps;++i){
        double u1= du(u,v,t);
        double v1= dv(u,v,t);
        double u2= du(u+0.5*u1*dt, v+0.5*v1*dt,t+0.5*dt);
        double v2= dv(u+0.5*u1*dt, v+0.5*v1*dt,t+0.5*dt);
        double u3= du(u+0.5*u2*dt, v+0.5*v2*dt,t+0.5*dt);
        double v3= dv(u+0.5*u2*dt, v+0.5*v2*dt,t+0.5*dt);
        double u4= du(u+u3*dt, v+v3*dt,t+dt);
        double v4= dv(u+u3*dt, v+v3*dt,t+dt);

        u+=dt/6.0 * (u1+2*u2+2*u3+u4);
        v+=dt/6.0 * (v1+2*v2+2*v3+v4);
        t+=dt;

        outFile << t << "\t" << u << "\t" << v << std::endl;
    }

    outFile.close();
}

duffing_oscillator::duffing_oscillator(double d, double a, double b, double g, double w): d(d), a(a), b(b), g(g), w(w){}