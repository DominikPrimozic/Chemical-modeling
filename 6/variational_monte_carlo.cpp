#include <variational_monte_carlo.h>

electron_in_potential::electron_in_potential(){
    std::random_device rd;
    gen = std::mt19937(rd());
}



double electron_in_potential::trial_wavefunction(double r, double a){
    return std::pow(a/M_PI,1.0/4.0)*std::exp(-a*0.5*(r*r));
}

double electron_in_potential::local_energy(double r, double a){ //for parabolic only
    return a*0.5 + (0.5-a*a*0.5)*r*r;
}

bool electron_in_potential::metropolis(double& r, double a, double delta){
    std::uniform_real_distribution<> dis(0,1);
    double r_=r + delta*(dis(gen)-0.5)*2;
    double psi=trial_wavefunction(r,a);
    double psi_=trial_wavefunction(r_,a);
    if (dis(gen)<(psi_*psi_) / (psi*psi)){
        r=r_;
        return true;
    }
    return false;
}

/*
double electron_in_potential::run(double a, int steps, double delta, int discard, int sample_interval){
    double r=0;
    expected_r=0;
    expected_r2=0;
    double energy=0;
    int sampled=0;
    for (int i=0;i<discard;i++){
        metropolis(r,a,delta);
    }
    for (int i=0;i<steps;i++){
        metropolis(r,a,delta);
        if (i%sample_interval==0){
            energy+=local_energy(r,a);
            sampled++;
            expected_r+=r;
            expected_r2+=r*r;
        }    
    }
    expected_r/=sampled;
    expected_r2/=sampled;
    return energy/sampled;
}
*/
double electron_in_potential::run(double a, int steps, double delta, int discard, int sample_interval, int tune){
    double energy = 0.0;
    expected_r = 0.0;
    expected_r2 = 0.0;

    int sampled = 0;

    #pragma omp parallel
    {
        std::mt19937 thread_gen(gen()); // unique RNG per thread
        double local_r = 0;
        double local_energy_sum = 0;
        double local_r_sum = 0;
        double local_r2_sum = 0;
        int local_sampled = 0;

        int taccept = 0, ttry = 0;

        for (int i = 0; i < discard; ++i) {
            if (metropolis(local_r, a, delta)) taccept++;
            ttry++;
            if (ttry==tune){
                double A = double(taccept) / ttry;
                if (A > 0.4) {
                    delta *=1.05;
                } else {
                    delta *=0.95;
                }
                taccept = ttry=0;
            }
        }

        for (int i = 0; i < steps; ++i) {
            metropolis(local_r, a, delta);
            if (i % sample_interval == 0) {
                double e = local_energy(local_r, a);
                local_energy_sum += e;
                local_r_sum += local_r;
                local_r2_sum += local_r * local_r;
                ++local_sampled;
            }
        }

        #pragma omp critical
        {
            energy += local_energy_sum;
            expected_r += local_r_sum;
            expected_r2 += local_r2_sum;
            sampled += local_sampled;
        }
    }

    expected_r /= sampled;
    expected_r2 /= sampled;
    return energy / sampled;
}

double electron_in_potential::optimize_alpha(double a_min, double a_max, double delta, double tol, int steps,int discard, int sample, int tune){
    double r=(std::sqrt(5)-1)/2.0;
    double x1 = a_max - r * (a_max - a_min);
    double x2 = a_min + r * (a_max - a_min);
    double f1 = run(x1, steps, delta,discard,sample,tune);
    double f2 = run(x2, steps, delta,discard,sample,tune);
    while (std::abs(a_max - a_min) > tol) {
        if (f1 < f2) {
            a_max = x2;
            x2 = x1;
            f2 = f1;
            x1 = a_max - r * (a_max - a_min);
            f1 = run(x1, steps, delta,discard,sample,tune);
        } else {
            a_min = x1;
            x1 = x2;
            f1 = f2;
            x2 = a_min + r * (a_max - a_min);
            f2 = run(x2, steps, delta,discard,sample,tune);
        }
    }

    double best_a = (a_min + a_max) / 2;
    std::cout << "Optimal a: " << best_a << "\nEnergy: " << run(best_a, steps, delta,discard,sample,tune)<< "\nExpected position: "<<expected_r<< "\nExpected position2: "<<expected_r2 << "\n";
    return best_a;
}