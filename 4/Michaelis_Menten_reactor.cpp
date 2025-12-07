#include <Michaelis_Menten_reactor.h>


MMdr::MMdr(double c0, double c_0, double c_h, double D, double v_max, double K_m, double h, double time):
    c0(c0), c_0(c_0),c_h(c_h), D(D), v_max(v_max), K_m(K_m),h(h), time(time){
    Pe=D*time/(h*h);
    Da=v_max*time/(c0*(K_m+c0));
}

double MMdr::MM_reaction(double z_i){
    return Da*z_i / (1+z_i);
}


void MMdr::solver( double dt, double dy){
    z.resize(static_cast<int>(1 / dy),1); //z=c/c0, starting is then 1
    std::vector<double> ztemp=z;
    double t=0;
    while (t<=1){
        for (int i=1;i<z.size()-1;i++){
            ztemp[i]= ( Pe*dt/(dy*dy) * (z[i+1]-2*z[i]+z[i-1])  +z[i] - dt* MM_reaction(z[i]));
        }
        t+=dt;
        z=ztemp;
    }
}

std::ostream& print_vector(std::ostream& os, const std::vector<double>& vec) {
    if (!vec.empty()) {
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            os << vec[i] << "\t"; 
        }
        os << vec.back();  
    }
    return os;
    return os;  
}

void MMdr::print_solver(std::string path, double dt, double dy){ //dt here is reduced time t/T
    std::cout << "Stability check value: " << Pe * dt / (dy * dy) << std::endl;
    if (Pe * dt / (dy * dy) > 0.5) {
        throw std::runtime_error("Stability condition violated");
    }

    z.resize(static_cast<int>(1/dy)+1,1); //z=c/c0, starting is then 1
    z[0]=c_0/c0;
    z[z.size()-1]=c_h/c0;
    std::vector<double> ztemp=z;
    double t=0;
    std::ofstream outFile(path);
    outFile << t << "\t";        
    print_vector(outFile, z);  
    outFile << std::endl;
    while (t<=1){
        for (int i=1;i<z.size()-1;i++){
            ztemp[i]= ( Pe*dt/(dy*dy) * (z[i+1]-2*z[i]+z[i-1])  +z[i] - dt* MM_reaction(z[i]));
        }
        t+=dt;
        z=ztemp;
        outFile << t << "\t";        
        print_vector(outFile, z);  
        outFile << std::endl;
    }
}