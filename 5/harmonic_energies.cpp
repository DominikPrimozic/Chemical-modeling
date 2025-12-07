#include <eigensolver.h>
#include <harmonic_oscillator.h>



int main(int argc, char* argv[]) {
    int n = 10;
    double alpha = 0;

    if (argc == 3) {
        n = std::stoi(argv[1]);
        alpha = std::stod(argv[2]);
    } else if (argc != 1) {
        std::cerr << "Usage: harmonic.exe [n alpha L N]" << std::endl;
        return 1;
    }

    quantum_harmonic a(n);
    a.solve_states(alpha);
    return 0;
}