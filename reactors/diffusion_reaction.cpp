#include <diffusion_reaction.h>


diffusion_reactor::diffusion_reactor(double D, double k, double a, double H, double c0) : D(D), k(k), rxn_order(a), H(H), c0(c0){
    outpath= "output/reaction/isolated/" + std::to_string(D) + "_"+ std::to_string(k) + "_" + std::to_string(a) + "_" + std::to_string(H) + "_" + std::to_string(c0) + ".txt";
}


double diffusion_reactor::du(double x, double u, double v){
    return v;
}


double diffusion_reactor::dv(double x, double u, double v){
    return k/D * std::pow(u,rxn_order);
}

double diffusion_reactor::RK4_solver(double dx, double s){
    double u=s; //this needs to be guessed
    double v=0;
    double x=0;

    while (x<=H){
        double k1=dx*du(x,u,v); //u is c and v is dc/dx
        double l1=dx*dv(x,u,v);
        double k2=dx*du(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double l2=dx*dv(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double k3=dx*du(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double l3=dx*dv(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double k4=dx*du(x+dx,u+k3,v+l3);
        double l4=dx*dv(x+dx,u+k3,v+l3);

        x+=dx;
        u+=1/6.0 * (k1+2*k2+2*k3+k4);
        v+=1/6.0 * (l1+2*l2+2*l3+l4);
    }
    
    return u;
}

//root function
double diffusion_reactor::root_function(double s,double dx){
    return RK4_solver(dx,s)-c0;
}

//root solver 
double diffusion_reactor::root_solver(double s, double h,double dx){
    //newton
    double s1;
    for (size_t i=0;i<10000;++i){
        double F = root_function(s,dx);
        double dF = (root_function(s+h,dx) - root_function(s-h,dx))/ (2*h);

        if (std::abs(dF) < 1e-10) {
            throw std::runtime_error("Derivative too small");
        }

        s1=s -  F/dF;
        if (std::abs(s-s1)<1e-5) return s1;
        s=s1;
    }
    throw std::runtime_error("Did not converge");
}

void diffusion_reactor::final_RK4_solver(double dx, double s){
    double u=s; //this was guessed
    double v=0; //nonporous wall
    double x=0;
    std::ofstream outFile(outpath);
    while (x<=H){
        outFile << x << "\t" << u<<"\t"<< v << std::endl;
        double k1=dx*du(x,u,v); //u is c and v is dc/dx
        double l1=dx*dv(x,u,v);
        double k2=dx*du(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double l2=dx*dv(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double k3=dx*du(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double l3=dx*dv(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double k4=dx*du(x+dx,u+k3,v+l3);
        double l4=dx*dv(x+dx,u+k3,v+l3);

        x+=dx;
        u+=1.0/6.0 * (k1+2*k2+2*k3+k4);
        v+=1.0/6.0 * (l1+2*l2+2*l3+l4);
    }
    outFile << x << "\t" << u<<"\t"<< v << std::endl;
    outFile.close();
}

void diffusion_reactor::solver(double dx, double s0, double h){
    double s=root_solver(s0,h,dx);
    final_RK4_solver(dx,s);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////temperature differnetials


temperature_diffusion_reactor::temperature_diffusion_reactor(double D, double k, double a, double H, double c0, double l, double rheat, double T0, double q, double Ea, double Ed) 
: D(D), k(k), rxn_order(a), H(H), c0(c0), heat(rheat), l(l), T0(T0), q(q), Ea(Ea), Ed(Ed){
    outpath= "output/reaction/flux/" + std::to_string(D) + "_"+  std::to_string(a) + "_" + std::to_string(c0) + "_"+
                                 std::to_string(T0) + "_"+ std::to_string(q) + ".txt";

    R=8.314;
}

double temperature_diffusion_reactor::reaction_constant(double w){
    return k*std::exp(-Ea/(R*w));
}

double temperature_diffusion_reactor::diffusion_constant(double w){
    return D*std::exp(-Ed/(R*w));
}


double temperature_diffusion_reactor::du(double x, double u, double v, double w, double z){
    return v;
}


double temperature_diffusion_reactor::dv(double x, double u, double v, double w, double z){
    return reaction_constant(w)/diffusion_constant(w) * std::pow(u,rxn_order);
}

double temperature_diffusion_reactor::dw(double x, double u, double v, double w, double z){
    return z;
}


double temperature_diffusion_reactor::dz(double x, double u, double v, double w, double z){
    return -heat/l * u;
}

double temperature_diffusion_reactor::RK4_solver_flux(double dx, double s){
    double u=s; //this needs to be guessed
    double v=0;
    double w=T0;
    double z=q;
    double x=0;

    while (x<=H){
        double k1=dx*du(x,u,v,w,z); //u is c and v is dc/dx, w is T and z is dT/dx
        double l1=dx*dv(x,u,v,w,z);
        double j1=dx*dw(x,u,v,w,z);
        double h1=dx*dz(x,u,v,w,z);
        double k2 = dx * du(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double l2 = dx * dv(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double j2 = dx * dw(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double h2 = dx * dz(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double k3 = dx * du(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double l3 = dx * dv(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double j3 = dx * dw(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double h3 = dx * dz(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double k4 = dx * du(x + dx, u + k3, v + l3, w + j3, z + h3);
        double l4 = dx * dv(x + dx, u + k3, v + l3, w + j3, z + h3);
        double j4 = dx * dw(x + dx, u + k3, v + l3, w + j3, z + h3);
        double h4 = dx * dz(x + dx, u + k3, v + l3, w + j3, z + h3);

        x+=dx;
        u+=1/6.0 * (k1+2*k2+2*k3+k4);
        v+=1/6.0 * (l1+2*l2+2*l3+l4);
        w+=1/6.0 * (j1+2*j2+2*j3+j4);
        z+=1/6.0 * (h1+2*h2+2*h3+h4);
    }
    
    return u;
}

//root function
double temperature_diffusion_reactor::root_function_flux(double s,double dx){
    return RK4_solver_flux(dx,s)-c0;
}

//root solver 
double temperature_diffusion_reactor::root_solver_flux(double s, double h,double dx){
    //newton
    double s1;
    for (size_t i=0;i<10000;++i){
        double F = root_function_flux(s,dx);
        double dF = (root_function_flux(s+h,dx) - root_function_flux(s-h,dx))/ (2*h);

        if (std::abs(dF) < 1e-10) {
            throw std::runtime_error("Derivative too small");
        }

        s1=s -  F/dF;
        if (std::abs(s-s1)<1e-5) return s1;
        s=s1;
    }
    throw std::runtime_error("Did not converge");
}

void temperature_diffusion_reactor::final_RK4_solver_flux(double dx, double s){
    double u=s; //this needs to be guessed
    double v=0;
    double w=T0;
    double z=q;
    double x=0;
    std::ofstream outFile(outpath);
    while (x<=H){
        outFile << x << "\t" << u<<"\t"<< v <<"\t"<<w<<"\t"<<z<< std::endl;
        double k1=dx*du(x,u,v,w,z); //u is c and v is dc/dx, w is T and z is dT/dx
        double l1=dx*dv(x,u,v,w,z);
        double j1=dx*dw(x,u,v,w,z);
        double h1=dx*dz(x,u,v,w,z);
        double k2 = dx * du(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double l2 = dx * dv(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double j2 = dx * dw(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double h2 = dx * dz(x + dx / 2.0, u + k1 / 2.0, v + l1 / 2.0, w + j1 / 2.0, z + h1 / 2.0);
        double k3 = dx * du(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double l3 = dx * dv(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double j3 = dx * dw(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double h3 = dx * dz(x + dx / 2.0, u + k2 / 2.0, v + l2 / 2.0, w + j2 / 2.0, z + h2 / 2.0);
        double k4 = dx * du(x + dx, u + k3, v + l3, w + j3, z + h3);
        double l4 = dx * dv(x + dx, u + k3, v + l3, w + j3, z + h3);
        double j4 = dx * dw(x + dx, u + k3, v + l3, w + j3, z + h3);
        double h4 = dx * dz(x + dx, u + k3, v + l3, w + j3, z + h3);

        x+=dx;
        u+=1/6.0 * (k1+2*k2+2*k3+k4);
        v+=1/6.0 * (l1+2*l2+2*l3+l4);
        w+=1/6.0 * (j1+2*j2+2*j3+j4);
        z+=1/6.0 * (h1+2*h2+2*h3+h4);
    }
    outFile << x << "\t" << u<<"\t"<< v <<"\t"<<w<<"\t"<<z<< std::endl;
    outFile.close();
}

void temperature_diffusion_reactor::solver_flux(double dx, double s0, double h){
    double s=root_solver_flux(s0,h,dx);
    final_RK4_solver_flux(dx,s);
}

//////////////////////////////////////////////////////////////////// temp profile

temperature_diffusion_reactor::temperature_diffusion_reactor(double D, double k, double a, double H, double c0, double l,double T0, double TH,double Ea, double Ed) 
    : D(D), k(k), rxn_order(a), H(H), c0(c0), l(l), T0(T0), TH(TH), Ea(Ea), Ed(Ed){

    outpath= "output/reaction/profile/" + std::to_string(D) + "_"+  std::to_string(a) + "_" + std::to_string(c0) + "_"+
                                 std::to_string(T0) + "_"+ std::to_string(TH) + ".txt";

    R=8.314;
}


double temperature_diffusion_reactor::temperature_profile(double x){
    return (TH-T0)/H * x + T0;
}

double temperature_diffusion_reactor::du(double x, double u, double v){
    return v;
}


double temperature_diffusion_reactor::dv(double x, double u, double v){
    return reaction_constant(temperature_profile(x))/diffusion_constant(temperature_profile(x)) * std::pow(u,rxn_order);
}

double temperature_diffusion_reactor::RK4_solver(double dx, double s){
    double u=s; //this needs to be guessed
    double v=0;
    double x=0;

    while (x<=H){
        double k1=dx*du(x,u,v); //u is c and v is dc/dx
        double l1=dx*dv(x,u,v);
        double k2=dx*du(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double l2=dx*dv(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double k3=dx*du(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double l3=dx*dv(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double k4=dx*du(x+dx,u+k3,v+l3);
        double l4=dx*dv(x+dx,u+k3,v+l3);

        x+=dx;
        u+=1/6.0 * (k1+2*k2+2*k3+k4);
        v+=1/6.0 * (l1+2*l2+2*l3+l4);
    }
    
    return u;
}

//root function
double temperature_diffusion_reactor::root_function(double s,double dx){
    return RK4_solver(dx,s)-c0;
}

//root solver 
double temperature_diffusion_reactor::root_solver(double s, double h,double dx){
    //newton
    double s1;
    for (size_t i=0;i<10000;++i){
        double F = root_function(s,dx);
        double dF = (root_function(s+h,dx) - root_function(s-h,dx))/ (2*h);

        if (std::abs(dF) < 1e-10) {
            throw std::runtime_error("Derivative too small");
        }

        s1=s -  F/dF;
        if (std::abs(s-s1)<1e-5) return s1;
        s=s1;
    }
    throw std::runtime_error("Did not converge");
}

void temperature_diffusion_reactor::final_RK4_solver(double dx, double s){
    double u=s; //this was guessed
    double v=0; //nonporous wall
    double x=0;
    std::ofstream outFile(outpath);
    while (x<=H){
        outFile << x << "\t" << u<<"\t"<< v << std::endl;
        double k1=dx*du(x,u,v); //u is c and v is dc/dx
        double l1=dx*dv(x,u,v);
        double k2=dx*du(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double l2=dx*dv(x+dx/2.0,u+k1/2.0,v+l1/2.0);
        double k3=dx*du(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double l3=dx*dv(x+dx/2.0,u+k2/2.0,v+l2/2.0);
        double k4=dx*du(x+dx,u+k3,v+l3);
        double l4=dx*dv(x+dx,u+k3,v+l3);

        x+=dx;
        u+=1.0/6.0 * (k1+2*k2+2*k3+k4);
        v+=1.0/6.0 * (l1+2*l2+2*l3+l4);
    }
    outFile << x << "\t" << u<<"\t"<< v << std::endl;
    outFile.close();
}

void temperature_diffusion_reactor::solver(double dx, double s0, double h){
    double s=root_solver(s0,h,dx);
    final_RK4_solver(dx,s);
}