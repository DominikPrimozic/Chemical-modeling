#include <eigensolver.h>
#include <harmonic_oscillator.h>



int main(int argc, char* argv[]) {
    int n = 10;
    double alpha = 0;
    double L = 5.0;
    int N = 1000;

    if (argc == 5) {
        n = std::stoi(argv[1]);
        alpha = std::stod(argv[2]);
        L = std::stod(argv[3]);
        N = std::stoi(argv[4]);
    } else if (argc != 1) {
        std::cerr << "Usage: harmonic.exe [n alpha L N]" << std::endl;
        return 1;
    }

    quantum_harmonic a(n);
    a.get_states(alpha, L, N);
    return 0;
}