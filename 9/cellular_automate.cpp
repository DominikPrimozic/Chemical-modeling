#include <cellular_automate.h>

grid_vector::grid_vector(int Nx_, int Ny_) : Nx(Nx_), Ny(Ny_) {
    g.resize(Nx * Ny, 0.0);
}
int& grid_vector::operator()(int x, int y) {
    return g[idx(x, y, Nx)];
}
int grid_vector::operator()(int x, int y) const {
    return g[idx(x, y, Nx)];
}


std::vector<int> BZ_reaction::count_excited_neighbors(int x, int y, int R){
    std::vector<int> count(3,0); //a,b,S
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0) continue;
            int nx=x+dx;
            int ny=y+dy;
            if (nx < 0 || nx >= current.Nx || ny < 0 || ny >= current.Ny) continue;

            if (current(nx, ny)>0 && current(nx, ny)<R) count[0]++;
            else if (current(nx, ny)==R) count[1]++;
            count[2]+=current(nx, ny);
        }
    }
    return count;
}
BZ_reaction::BZ_reaction(int Nx, int Ny) : current(Nx, Ny), next(Nx, Ny){};


void BZ_reaction::ca_step(int R, int k1, int k2, int g) {
    #pragma omp parallel for
    for (int y = 0; y < current.Ny; ++y) {
        for (int x = 0; x < current.Nx; ++x) {
            int state = current(x, y);
            std::vector<int> a_b_S = count_excited_neighbors(x,y,R);
            if (state == 0) {
                next(x,y)= (int)(a_b_S[0]/k1) + (int)(a_b_S[1]/k2);
            } else if (state>=R) {
                next(x, y) = 0;
            } else {
                next(x, y) = (int)(a_b_S[2]/(a_b_S[0]+1)) + g; 
            }
        }
    }
    std::swap(current, next);
}

void BZ_reaction::write_frame(int frame_number, int R) {
    std::ostringstream filename;
    filename << "CA/frame_" << std::setw(4) << std::setfill('0') << frame_number << ".ppm";
    std::ofstream file(filename.str(), std::ios::binary);

    file << "P3\n" << current.Nx << " " << current.Ny << "\n255\n";

    for (int y = 0; y < current.Ny; ++y) {
        for (int x = 0; x < current.Nx; ++x) {
            int s = current(x, y);
            int r = 0, g = 0, b = 0;

            if (s == 0) {
                r = g = b = 0;
            } else if (s == 1) {
                r = 255; g = 0; b = 0;
            } else {
                float t = float(s - 2) / float(R - 2);
                t = std::max(0.0f, std::min(1.0f, t));  // Clamp t between 0 and 1
                r = int((1 - t) * 255);
                g = int(t * 255);
                b = 128;
            }

            file << r << " " << g << " " << b << " ";
        }
        file << "\n";
    }
}


void BZ_reaction::run(int steps, int R, int k1, int k2, int g) {
        for (int t = 0; t < steps; ++t) {
            ca_step(R,k1,k2,g);
            write_frame(t, R);

            if (t==steps/2){
                int active = 0;
                for (int y = 0; y < current.Ny; ++y)
                    for (int x = 0; x < current.Nx; ++x)
                        if ((int)current(x, y) > 0) active++;

                std::cout << "Active cells at step " << t << ": " << active << "\n";
                }
        }
    }

void BZ_reaction::set_state(int x, int y, int value) {
        current(x, y) = value;
    }

void BZ_reaction::random_initialize(double p_excited, double p_refractory, int R) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int y = 0; y < current.Ny; ++y) {
        for (int x = 0; x < current.Nx; ++x) {
            double r = dis(gen);
            if (r < p_excited) {
                current(x, y) = 1;
            } else if (r < p_excited + p_refractory) {
                current(x, y) = 2 + int(dis(gen) * (R - 2));
            } else {
                current(x, y) = 0;
            }
        }
    }
}