#include <lbm.h>
#include <filesystem> 

namespace fs = std::filesystem;

int main() {
    // -------------------------
    // Simulation grid parameters
    // -------------------------
    int Nx = 100;            // Lattice width
    int Ny = 50;             // Lattice height (channel height = H)
    double H = Ny;           // Channel height in lattice units

    // -------------------------
    // Flow configuration
    // -------------------------
    double Re = 1000.0;       // Desired Reynolds number
    double Umax = 0.05;      // Target maximum velocity (for stability: Umax << 1)

    // -------------------------
    // LBM unit calculations
    // -------------------------
    double nu = (Umax * H) / Re;       // Kinematic viscosity from Re = U L / ν
    double tau = 3.0 * nu + 0.5;       // LBM relaxation time
    double gravity = (8.0 * nu * Umax) / (H * H);  // Gravity to match Poiseuille max velocity

    // Print computed parameters
    std::cout << "Computed LBM Parameters:\n";
    std::cout << " - nu   = " << nu << "\n";
    std::cout << " - tau  = " << tau << "\n";
    std::cout << " - g    = " << gravity << "\n";
    std::cout << " - H    = " << H << "\n";
    std::cout << " - Re   = " << Re << "\n";

    // -------------------------
    // Create output folder if not exists
    // -------------------------
     fs::path outdir("output");
    if (fs::exists(outdir)) {
        for (const auto& entry : fs::directory_iterator(outdir)) {
            fs::remove(entry.path());  // Delete each file
        }
    } else {
        fs::create_directory(outdir);  // Create if missing
    }

    // -------------------------
    // Initialize simulation with parameters
    // -------------------------
    D2Q9_2D sim(Nx, Ny,Re,Umax,H);

    // -------------------------
    // Run simulation
    // -------------------------
    sim.run();

    return 0;
}