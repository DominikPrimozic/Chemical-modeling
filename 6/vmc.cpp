#include <variational_monte_carlo.h>
#include <cstdlib>

int main(int argc, char* argv[]) {
    double a_min=0.1, a_max=2.0;
    double delta=1;
    int steps=100000, discard=10000, sample=1000, tune=100;
    if (argc > 1) a_min    = std::atof(argv[1]);
    if (argc > 2) a_max    = std::atof(argv[2]);
    if (argc > 3) delta    = std::atof(argv[3]);
    if (argc > 4) steps    = static_cast<int>(std::atof(argv[4]));
    if (argc > 5) discard  = static_cast<int>(std::atof(argv[5]));
    if (argc > 6) sample   = static_cast<int>(std::atof(argv[6]));
    if (argc > 7) tune     = static_cast<int>(std::atof(argv[7]));

    if (argc > 8) {
        std::cerr << "Usage: " << argv[0]
                  << " [a_min a_max delta steps discard sample tune]\n";
        return 1;
    }

    electron_in_potential system;
    system.optimize_alpha(a_min,a_max,delta,1e-5,steps,discard,sample,tune);
    return 0;
}