#include <molecular_dynamics.h>



LJ_MD::LJ_MD(int N, double density, double placement_radius, std::string pat, bool log, int dim) : N(N), dim(dim), log(log), pat(pat){
    random_box placement(N,density,dim,placement_radius);
    r=placement.getItems();
    L=placement.getBox();
    v.resize(N,std::vector<double>(dim));
    f.resize(N,std::vector<double>(dim));
    force_temp.resize(dim);
    f_temp.resize(N,std::vector<double>(dim));
    difr.resize(dim);
}

void velocity_sampling(double T, std::vector<std::vector<double>>& v){
    //technically we need *m/kb in the exponent
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, std::sqrt(T));
    int dim=v[0].size();
    int N=v.size();
    std::vector<double> v_c(dim,0);
    for (int p=0;p<N;p++){
        for (int d=0;d<dim;d++){
            double temp=dist(gen);
            v[p][d]=temp;
            v_c[d]+=temp;

        }
    }
    //keeping center of mass
    for (int d=0;d<dim;d++){
        v_c[d]/=(double)N;
    }  
    //shift
    #pragma omp parallel for
    for (int p=0;p<N;p++){
        for (int d=0;d<dim;d++){
            v[p][d]-=v_c[d];
        }
    }
}

double LJ_MD::distance2(std::vector<double>& dr, const std::vector<double>& p1, const std::vector<double>& p2){
    double r=0;
    for (size_t i=0;i<dim;i++){
        double delta =(p1[i]-p2[i]);
        delta -= L * std::round(delta / L);
        dr[i] = delta;
        r+=delta*delta;
    }
    return r;  
}

double LJ_MD::pair_force(double r2, double sigma, double epsilon){
    if (r2 < 1e-12) r2 = 1e-12;
    double rm2=sigma*sigma/r2;
    double rm6=rm2*rm2*rm2;
    
    
    return (24*epsilon*(2*rm6*rm6 - rm6) / r2); //normalization
}

void LJ_MD::update_force(){
    set_zero(f); //so we do not reallocate
    for (int i=0;i<N;i++){
        for (int j=i+1;j<N;j++){
            double r2=distance2(difr,r[i],r[j]);
            if (r2<r_cut2){
                force_temp=(pair_force(r2))*(difr);
                f[i]+=force_temp;
                f[j]-=force_temp;
            }
        }
    }
}

void LJ_MD::velocity_verlet(double dt){
    double m=1.0;
    #pragma omp parallel for
    for (int i=0;i<N;i++){
        for (int d = 0; d < dim; d++) {
            r[i][d] += v[i][d]*dt + (0.5/m) * f[i][d]*dt*dt;
            r[i][d] -= L * std::round(r[i][d] / L);
        }
        f_temp[i]= f[i];
    }
    
    if (track_positions) write_positions_binary(position);
    update_force();
    #pragma omp parallel for
    for (int i=0;i<N;i++){
        v[i]+=0.5*dt/m * (f_temp[i] + f[i]);
    }
    if (track_velocities) write_velocities_binary(velocity);
}

double LJ_MD::kinetic_energy() {
    double m=1;
    double kinetic_energy=0; 
    #pragma omp parallel for reduction(+:kinetic_energy)
    for (int p=0;p<N;p++){
        double v2=0;
        for (int d=0;d<dim;d++){
            v2+=v[p][d]*v[p][d];
        }
        kinetic_energy+=v2;
    }
    return kinetic_energy*0.5*m;
}

double LJ_MD::measure_temperature(double KE) {
    return (2.0 * KE) / (dim * N); 
}

double LJ_MD::pressure() {
    double volume = std::pow(L, dim);
    double virial_sum = 0.0;
    double r2;
    
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            r2 = distance2(difr, r[i], r[j]);
            if (r2 < r_cut2) {
                double f_scalar = pair_force(r2);  
                for (int d = 0; d < dim; ++d) {
                    virial_sum += difr[d] * f_scalar * difr[d];
                }
            }
        }
    }

    double T_inst = measure_temperature(kinetic_energy());
    double ideal = N * T_inst / volume;
    double virial = virial_sum / (dim * volume);

    return ideal + virial;
}



double LJ_MD::potential_energy(double epsilon, double sigma){ 
    double potential=0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double r2 = distance2(difr, r[i], r[j]);
            if (r2 < r_cut2) {
                double rm2 = sigma * sigma / r2; 
                double rm6 = rm2 * rm2 * rm2;
                if (r2 < 1e-12) {
                    r2= 1e-12;  
                }
                potential+=4*epsilon*(rm6*rm6-rm6);
            }
        }
    }
    return potential;
}

double LJ_MD::total_energy() {
    return kinetic_energy() + potential_energy();
}

void LJ_MD::rescale_thermostat(){
     
    double T_actual=measure_temperature(kinetic_energy());

    double lambda=std::sqrt(T/T_actual);
    #pragma omp parallel for
    for (int p=0;p<N;p++){
        for (int d=0;d<dim;d++){
            v[p][d]*=lambda;
        }
    }
}

void LJ_MD::run(double setT, int equilibration_steps, int steps, double dt, int thermostat_interval, int sampling_interval, double rc){
    fs::create_directories(pat);
    fs::create_directories(pat + "pos/");
    fs::create_directories(pat + "vel/");
    
    T=setT;
    r_cut2=rc*rc;
    write_metadata(pat + "metadata.txt", equilibration_steps, steps,dt);
    velocity_sampling(T,v);
    update_force();
    for (int step = 0; step < equilibration_steps; step++){
        velocity_verlet(dt);
        /*
        if (step % 1000 == 0) {
            // 1) kinetic
            double KE = kinetic_energy();
            double T_inst = measure_temperature(KE);
    
            std::cout << "EQ step=" << step
                      << "  T=" << T_inst
                      << "\n";
        }*/
        if (step  % thermostat_interval == 0) {
            double KE = kinetic_energy();
            double T_inst = measure_temperature(KE);
    
            std::cout << "EQ step=" << step
                      << "  T=" << T_inst
                      << "\n";

            rescale_thermostat();
        }
    }
    std::ofstream thermofile( pat + "thermo.txt");
    thermofile << "# Time\tTemperature\tKinetic_E\tPotential_E\tTotal_E\tPressure\n";
    for (int step = 0; step < steps; step++) {
        if (log && step % sampling_interval == 0) {
            track_positions = true;
            track_velocities = true;
    
            position = pat + "pos/" + std::to_string(step * dt) + ".bin";
            velocity = pat + "vel/" + std::to_string(step * dt) + ".bin";
        }
    
        velocity_verlet(dt);
    
        if (step % sampling_interval == 0) {
            double kinE = kinetic_energy();
            double potE = potential_energy();
            double totE = kinE + potE;
            double temp = measure_temperature(kinE);
            double pres = pressure();
    
            thermofile << step * dt << "\t"
                       << temp << "\t"
                       << kinE << "\t"
                       << potE << "\t"
                       << totE << "\t"
                       << pres << "\n";
    
            rescale_thermostat();
    
            track_positions = false;
            track_velocities = false;
        }
    }
    thermofile.close();
}


void LJ_MD::write_positions_binary(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    for (const auto& p : r) {
        out.write(reinterpret_cast<const char*>(p.data()), sizeof(double) * p.size());
    }
    out.close();
}

void LJ_MD::write_velocities_binary(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    for (const auto& p : v) {
        out.write(reinterpret_cast<const char*>(p.data()), sizeof(double) * p.size());
    }
    out.close();
}

void LJ_MD::write_metadata(const std::string& filename, int equilibration_steps,int steps, double dt) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Failed to open metadata file: " << filename << std::endl;
        return;
    }

    out << "# LJ MD Simulation Metadata\n";
    out << "N = " << N << "\n";
    out << "dim = " << dim << "\n";
    out << "L = " << L << "\n";
    out << "timestep = " << dt << "\n";
    out << "equilibration steps = " << equilibration_steps << "\n";
    out << "simulation steps = " << steps << "\n";
    out << "temperature = " << T << "\n";
    out << "r_cut = " << std::sqrt(r_cut2) << "\n";
    out.close();
}