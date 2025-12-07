#include <cluster_dynamics.h>



LJ_cluster::LJ_cluster(std::string pat, bool log, int dim) : dim(dim), log(log), pat(pat){
    //9 particles of radius 1
    double spacing=2*1;
    if (dim==2){
        r.push_back({0, 0});
        r.push_back({-spacing, -spacing});
        r.push_back({spacing, -spacing});
        r.push_back({-spacing * 2, -2 * spacing});
        r.push_back({0, -2 * spacing});
        r.push_back({spacing * 2, -2 * spacing});
        r.push_back({-spacing, -3 * spacing});
        r.push_back({spacing, -3 * spacing});
        r.push_back({0, -4 * spacing});
    }
    if (dim==3){
        r.push_back({-spacing, -spacing, -spacing}); // Particle 1 
        r.push_back({spacing, -spacing, -spacing});  // Particle 2 
        r.push_back({-spacing, spacing, -spacing});  // Particle 3 
        r.push_back({spacing, spacing, -spacing});   // Particle 4 
        r.push_back({-spacing * 2, -2 * spacing, 0});  // Particle 5 
        r.push_back({0, -2 * spacing, 0});             // Particle 6
        r.push_back({spacing * 2, -2 * spacing, 0});   // Particle 7 
        r.push_back({-spacing, -spacing, 0});          // Particle 8 
        r.push_back({0, -spacing, 0});                 // Particle 9 
        r.push_back({spacing, -spacing, 0});           // Particle 10 
        r.push_back({-spacing * 2, 0, 0});             // Particle 11 
        r.push_back({0, 0, 0});                       // Particle 12 
        r.push_back({spacing * 2, 0, 0});             // Particle 13 
        r.push_back({-spacing, -spacing, spacing});   // Particle 14 
        r.push_back({spacing, -spacing, spacing});    // Particle 15 
        r.push_back({-spacing, spacing, spacing});    // Particle 16 
        r.push_back({spacing, spacing, spacing});     // Particle 17 
    }
    N=r.size();
    v.resize(N,std::vector<double>(dim));
    f.resize(N,std::vector<double>(dim));
    force_temp.resize(dim);
    f_temp.resize(N,std::vector<double>(dim));
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

double LJ_cluster::distance2(const std::vector<double>& p1, const std::vector<double>& p2){
    double r=0;
    for (size_t i=0;i<dim;i++){
        double delta =(p1[i]-p2[i]);
        r+=delta*delta;
    }
    return r;  
}

double LJ_cluster::pair_force(double r2, double sigma, double epsilon){
    if (r2 < 1e-12) r2 = 1e-12;
    double rm2=sigma*sigma/r2;
    double rm6=rm2*rm2*rm2;
    
    
    return (24*epsilon*(2*rm6*rm6 - rm6) / r2); //normalization
}

void LJ_cluster::update_force(){
    set_zero(f); //so we do not reallocate
    for (int i=0;i<N;i++){
        for (int j=i+1;j<N;j++){
            double r2=distance2(r[i],r[j]);
            force_temp=(pair_force(r2))*(r[i]-r[j]);
            f[i]+=force_temp;
            f[j]-=force_temp;
            
        }
    }
}

void LJ_cluster::velocity_verlet(double dt){
    double m=1.0;
    #pragma omp parallel for
    for (int i=0;i<N;i++){
        for (int d = 0; d < dim; d++) {
            r[i][d] += v[i][d]*dt + (0.5/m) * f[i][d]*dt*dt;
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

double LJ_cluster::kinetic_energy() {
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

double LJ_cluster::measure_temperature(double KE) {
    return (2.0 * KE) / (dim * N); 
}


double LJ_cluster::potential_energy(double epsilon, double sigma){ 
    double potential=0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double r2 = distance2(r[i], r[j]);
            double rm2 = sigma * sigma / r2; 
            double rm6 = rm2 * rm2 * rm2;
            if (r2 < 1e-12) {
                r2= 1e-12;  
            }
            potential+=4*epsilon*(rm6*rm6-rm6);
            
        }
    }
    return potential;
}

double LJ_cluster::total_energy() {
    return kinetic_energy() + potential_energy();
}

void LJ_cluster::rescale_thermostat(){
     
    double T_actual=measure_temperature(kinetic_energy());

    double lambda=std::sqrt(T/T_actual);
    #pragma omp parallel for
    for (int p=0;p<N;p++){
        for (int d=0;d<dim;d++){
            v[p][d]*=lambda;
        }
    }
}

void LJ_cluster::run(double setT, int equilibration_steps, int steps, double dt, int thermostat_interval, int sampling_interval){
    fs::create_directories(pat);
    fs::create_directories(pat + "pos/");
    fs::create_directories(pat + "vel/");
    
    T=setT;
    write_metadata(pat + "metadata.txt", equilibration_steps, steps,dt);
    //velocity_sampling(T,v);
    update_force();
    for (int step = 0; step < equilibration_steps; step++){
        velocity_verlet(dt);
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
    
            thermofile << step * dt << "\t"
                       << temp << "\t"
                       << kinE << "\t"
                       << potE << "\t"
                       << totE << "\n";
    
            rescale_thermostat();
    
            track_positions = false;
            track_velocities = false;
        }
    }
    thermofile.close();
}


void LJ_cluster::write_positions_binary(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    for (const auto& p : r) {
        out.write(reinterpret_cast<const char*>(p.data()), sizeof(double) * p.size());
    }
    out.close();
}

void LJ_cluster::write_velocities_binary(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    for (const auto& p : v) {
        out.write(reinterpret_cast<const char*>(p.data()), sizeof(double) * p.size());
    }
    out.close();
}

void LJ_cluster::write_metadata(const std::string& filename, int equilibration_steps,int steps, double dt) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Failed to open metadata file: " << filename << std::endl;
        return;
    }

    out << "# LJ MD Simulation Metadata\n";
    out << "N = " << N << "\n";
    out << "dim = " << dim << "\n";
    out << "timestep = " << dt << "\n";
    out << "equilibration steps = " << equilibration_steps << "\n";
    out << "simulation steps = " << steps << "\n";
    out << "temperature = " << T << "\n";
    out.close();
}

void LJ_cluster::random_initialize(double box_size) {
    r.clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-box_size, box_size);
    for (int i = 0; i < N; i++) {
        std::vector<double> pos(dim);
        for (int d = 0; d < dim; d++) {
            pos[d] = dis(gen);
        }
        r.push_back(pos);
    }
}

void LJ_cluster::steepest_descent_minimize(double dt, int max_steps) {
    update_force();
    for (int step = 0; step < max_steps; ++step) {
        double max_force = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int d = 0; d < dim; ++d) {
                r[i][d] += dt * f[i][d]; // small step in force direction
                max_force = std::max(max_force, std::abs(f[i][d]));
            }
        }
        update_force();
        if (max_force < 1e-5) { // convergence criterion
            break;
        }
    }
}