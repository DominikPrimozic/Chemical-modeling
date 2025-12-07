#include <monte_carlo.h>



Lennard_Jones::Lennard_Jones(int N, double density, double T, double e, double s, int dim, double rc) : N(N),density(density),T(T),e(e),s(s), dimensions(dim){
    random_box initial(N,density,dimensions,s);
    particles=initial.getItems();
    L=initial.getBox();
    trial_position.resize(dimensions);

    r_cutoff=rc;
    r_cutoff2=rc*rc;
    s2=s*s;

    std::random_device rd;
    gen = std::mt19937(rd());

    
}


double Lennard_Jones::distance2(std::vector<double>& p1, std::vector<double>& p2){
    double r=0;
    for (size_t i=0;i<dimensions;i++){
        double delta =(p1[i]-p2[i]);
        //delta = delta - std::floor((delta + L/2.0) / L) * L;
        delta -= L * std::round(delta / L);
        r+=delta*delta;
    }
    return r;  
}

double Lennard_Jones::LJ_pair(double r2){
    const double epsilon = 1e-12;
    if (r2 < epsilon) {
        return 1e12;  // huge repulsion for overlapping particles
    }
    return 4*e*(std::pow(s2/r2,6.0)-std::pow(s2/r2,3.0));
}



double Lennard_Jones::LJ_potential(int i,std::vector<double>& i_coords){ //i<N is a must
    double potential=0;
    for (size_t j = 0; j < N; ++j) {
        if (j == i) continue; 
        double r2=distance2(i_coords, particles[j]);
        if (r2 < r_cutoff*r_cutoff){
            potential += LJ_pair(r2);
        }
    }
    return potential;
}

bool Lennard_Jones::metropolis(double dE){
    std::uniform_real_distribution<> mep(0,1);
    if (dE < 0) return true;  
    return mep(gen) < std::exp(-dE / T);
}


bool Lennard_Jones::move(double step_size){
    std::uniform_real_distribution<> dis(-step_size,step_size);
    std::uniform_int_distribution<> particle_picker(0,N-1);
    int particle_move=particle_picker(gen);
    trial_position = particles[particle_move];
    double E_i=LJ_potential(particle_move,trial_position);
    for (int i=0;i<dimensions;i++){
        trial_position [i]+=dis(gen);
        trial_position [i] -= L * std::round(trial_position [i] / L);
    }
    double E_new=LJ_potential(particle_move,trial_position);
    double dE=E_new-E_i;

    if (metropolis(dE)){
        particles[particle_move]=trial_position ;
        return true;
    }
    return false;
}



void Lennard_Jones::run_simulation(simulation_parameters& sp){
    if (sp.keep_log) {
        std::ofstream(sp.log_path + "step_size.txt", std::ios::trunc).close();
        std::ofstream(sp.log_path + "acceptance_log.txt", std::ios::trunc).close();
        std::ofstream(sp.energies + "equilibration_energy.txt", std::ios::trunc).close();
        std::ofstream(sp.energies + "sampled_energy.txt", std::ios::trunc).close();
    }
    std::filesystem::create_directories(sp.configurations);
    std::filesystem::create_directories(sp.energies);
    std::filesystem::create_directories(sp.log_path);
    E=compute_energy(); //starting energy
    equilibrate(sp);
    sampling(sp);

}



double Lennard_Jones::compute_energy(){
    double E_new=0;
    #pragma omp parallel for reduction(+:E_new) schedule(dynamic)
    for (int i=0;i<N;i++){
        for (int j=i+1;j<N;j++){
            double r2=distance2(particles[i], particles[j]);
            if (r2 < r_cutoff2){
                E_new += LJ_pair(r2);
            }
        }
    }
    return E_new;
}













void Lennard_Jones::equilibrate(simulation_parameters& sp) {
    sp.step_size = 0.1*s;
    double& ss = sp.step_size;
    int accepted = 0, attempted = 0, step = 0;
    double E_new = E;
    std::ofstream equil_energy_out(sp.energies + "equilibration_energy.txt");

    while (true) {
        bool success = move(ss);
        attempted++;
        if (success) accepted++;
        step++;

        if (step % sp.equil_steps_check == 0) {
            E_new = compute_energy();
            equil_energy_out << step << " " << E_new << "\n";

            if (step == sp.equil) {
                std::cout << "System equilibrated at step " << std::endl;
                std::cout << "Box length " <<L <<std::endl;
                std::cout << "Step length " <<ss <<std::endl;
                break;
            }
            E = E_new;
        }

        if (step % sp.tune == 0 && step > 0) {
            double accept_ratio = (double)accepted / attempted;
            if (accept_ratio > 0.6) ss *= 1.1;
            else if (accept_ratio < 0.3) ss *= 0.9;
            //std::cout << "Step length " <<ss <<std::endl;
            if (ss > L) ss = L;
            if (sp.keep_log) {
                std::ofstream step_log(sp.log_path + "step_size.txt", std::ios::app);
                step_log << step << " " << ss << "\n";
            }
            accepted = 0;
            attempted = 0;
        }
    }

    equil_energy_out.close();
}

void Lennard_Jones::sampling(simulation_parameters& sp) {
    double& ss = sp.step_size;
    int accepted = 0, attempted = 0;
    std::ofstream energy_out(sp.energies + "sampled_energy.txt");

    for (int step = 0; step < sp.steps; ++step) {
        bool success = move(ss);
        attempted++;
        if (success) accepted++;

        if (step % sp.sample_interval == 0) {
            sample_step(step, energy_out, sp);
            if (sp.keep_log) {
                std::ofstream log(sp.log_path + "acceptance_log.txt", std::ios::app);
                log << step << " " << sp.step_size << " " << (double)accepted / attempted  << "\n";
                //std::cout << "[Sampling] Step " << step
                //          << " | Acceptance: " << (double)accepted / attempted << "\n";
            }
        }
        
    }

    energy_out.close();
}

void Lennard_Jones::sample_step(int step, std::ofstream& out_energy, const simulation_parameters& sp) {
    double current_energy = compute_energy();
    out_energy << step << " " << current_energy << "\n";

    if (sp.save_configs) {
        std::string fname = sp.configurations + "config_step_" + std::to_string(step) + ".bin";
        write_particles_binary(fname);
    }
}

void Lennard_Jones::write_particles_binary(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    for (const auto& p : particles) {
        out.write(reinterpret_cast<const char*>(p.data()), sizeof(double) * p.size());
    }
    out.close();
}
