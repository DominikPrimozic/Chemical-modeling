#include <distributions.h>

boxmuller_double::boxmuller_double(double u1, double u2) {
    z1 = std::sqrt(-2 * std::log(u1)) * std::cos(2 * PI * u2);
    z2 = std::sqrt(-2 * std::log(u1)) * std::sin(2 * PI * u2);
}

int value_to_bin(double val, int bins, double min_val,double max_val){
    double bin_width=(max_val-min_val)/bins;
    int bin_index=static_cast<int>((val-min_val)/bin_width);
    if (bin_index < 0) {bin_index = 0;}
    if (bin_index >= bins) {bin_index = bins - 1;}

    return bin_index;
}

maxwell_distribution::maxwell_distribution(int N, std::string path, double T, double m,double min_val,double max_val) : T(T), m(m){
    k=1.38e-23;
    sigma=std::sqrt(k*T/m);
    generate_distribution(N,min_val,max_val,path);
    

}

void maxwell_distribution::generate_distribution(int N, double min_val,double max_val,std::string path){
    //int bins = static_cast<int>(1 + std::log2(N)); //shit
    int bins=500;

    x_bins.resize(bins, 0);
    v_bins.resize(bins, 0);
    E_bins.resize(bins, 0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0,1);
    double u1,u2,u3,u4;

    
    for (size_t i=0;i<N;++i){
        u1=dis(gen);
        u2=dis(gen);
        u3=dis(gen);
        u4=dis(gen);
        boxmuller_double boxmuller_double1(u1,u2);
        boxmuller_double boxmuller_double2(u3,u4);

        //for vx
        int current_bin_vx= value_to_bin( boxmuller_double1.z1*sigma,bins,min_val,max_val);
        x_bins[current_bin_vx]++;

        //for v
        double v=std::sqrt(boxmuller_double1.z1*boxmuller_double1.z1*sigma*sigma + boxmuller_double1.z2*boxmuller_double1.z2*sigma*sigma + boxmuller_double2.z1*boxmuller_double2.z1*sigma*sigma);
        int current_bin_v= value_to_bin( v,bins,0,max_val);
        v_bins[current_bin_v]++;

        //for E
        double E=0.5*m*v*v/1.6e-19;
        int current_bin_E= value_to_bin( E ,bins,0,m*0.5*max_val*max_val/1.6e-19);
        E_bins[current_bin_E]++;
    }

    std::string vx_path = path + "vx.txt";
    std::string v_path = path + "v.txt";
    std::string E_path = path + "E.txt";
    print_to_file(vx_path, bins, x_bins);
    print_to_file(v_path,  bins, v_bins);
    print_to_file(E_path,  bins, E_bins);
    
}

void maxwell_distribution::print_to_file(std::string path,int bins, std::vector<int>& bin_vector){
    std::ofstream file(path);

    for (int i = 0; i < bins; ++i) {
        file << i << "\t" << bin_vector[i] << std::endl;
    }
    file.close();
}