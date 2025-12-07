#include <lbm.h>

lattice_velocities::lattice_velocities(int Nx_, int Ny_) : Nx(Nx_), Ny(Ny_) {
    for (int i = 0; i < 9; ++i)
        f[i].resize(Nx * Ny, 0.0);
}

double& lattice_velocities::operator()(int i, int x, int y) {
    return f[i][idx(x, y, Nx)];
}
double lattice_velocities::operator()(int i, int x, int y) const {
    return f[i][idx(x, y, Nx)];
}
double& lattice_velocities::operator()(int i, int idx) {
    return f[i][idx];
}
double lattice_velocities::operator()(int i, int idx) const {
    return f[i][idx];
}

grid_vector::grid_vector(int Nx_, int Ny_) : Nx(Nx_), Ny(Ny_) {
    g.resize(Nx * Ny, 0.0);
}
double& grid_vector::operator()(int x, int y) {
    return g[idx(x, y, Nx)];
}
double grid_vector::operator()(int x, int y) const {
    return g[idx(x, y, Nx)];
}

void write_field_to_file(const grid_vector& field, const std::string& filename, int step) {
    std::ofstream out(filename + "_" + std::to_string(step) + ".dat");
    if (!out) {
        std::cerr << "Error: Could not open file " << filename << "_" << step << ".dat" << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(6);
    for (int y = 0; y < field.Ny; ++y) {
        for (int x = 0; x < field.Nx; ++x) {
            out << field(x, y) << " ";
        }
        out << "\n";
    }
}

void print_grid_vector(const grid_vector& rho) {
    for (int y = 0; y < rho.Ny; ++y) {
        for (int x = 0; x < rho.Nx; ++x) {
            std::cout << rho(x, y) << " ";
        }
        std::cout << std::endl;
    }
}

D2Q9_2D::D2Q9_2D(int Nx_, int Ny_, double Re, double U, double L) :
    f(Nx_, Ny_), f0(Nx_, Ny_), rho(Nx_, Ny_), vx(Nx_, Ny_), vy(Nx_, Ny_), obstacle_mask(Nx_, Ny_),
    gravity(1e-6), dt(1.0), Nx(Nx_), Ny(Ny_), dx(1.0)
{
    gravity = (8.0 * nu * U) / (L*L);
    nu = U * L / Re;
    tau = 3.0 * nu + 0.5;
    w = {4/9.0, 1/9.0, 1/9.0, 1/9.0, 1/9.0, 1/36.0, 1/36.0, 1/36.0, 1/36.0};
    cx = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    cy = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    opp = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            rho(x, y) = 1.0;
            vx(x, y) = U;   
            vy(x, y) = 0.0;
        }
    }
    // Initialize obstacle
    int obstacle_width = 5;
    int half_width = obstacle_width / 2;
    int start=Nx/3; //Nx/2
    int x_start = start - half_width;
    int x_end = start + half_width;

    for (int x = x_start; x <= x_end; ++x) {
        //for (int y = Ny / 4; y < 3 * Ny / 4; ++y) {
        for (int y = Ny / 2 - Ny / 8; y < Ny / 2 + Ny / 8; ++y){
            obstacle_mask(x, y) = 1.0;
        }
    }

    equilibrium_f();
    for (int i = 0; i < 9; ++i)
        for (int idx_ = 0; idx_ < Nx * Ny; ++idx_)
            f(i, idx_) = f0(i, idx_);
}

void D2Q9_2D::equilibrium_f() {
    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            double u_sq = vx(x, y) * vx(x, y) + vy(x, y) * vy(x, y);
            for (int i = 0; i < 9; ++i) {
                double cu = cx[i] * vx(x, y) + cy[i] * vy(x, y);
                f0(i, x, y) = w[i] * rho(x, y) * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * u_sq);
            }
        }
    }
}

void D2Q9_2D::collision() {
    for (int i = 0; i < 9; ++i) {
        for (int idx = 0; idx < Nx * Ny; ++idx) {
            f(i, idx) -= (1.0 / tau) * (f(i, idx) - f0(i, idx));
        }
    }
}

void D2Q9_2D::stream() {
    lattice_velocities f_new(Nx, Ny);
    for (int i = 0; i < 9; ++i) {
        int dx = cx[i];
        int dy = cy[i];
        for (int x = 0; x < Nx; ++x) {
            for (int y = 0; y < Ny; ++y) {
                int x_pull = (x - dx + Nx) % Nx;
                int y_pull = (y - dy + Ny) % Ny;

                f_new(i, x, y) = f(i, x_pull, y_pull);
            }
        }
    }
    std::swap(f, f_new);
}

void D2Q9_2D::bounce_boundaries() {
    for (int x = 0; x < Nx; ++x) {
        int y_top = Ny - 1;
        int y_bot = 0;

        f(4, x, y_top) = f0(2, x, y_top);
        f(7, x, y_top) = f0(5, x, y_top);
        f(8, x, y_top) = f0(6, x, y_top);

        f(2, x, y_bot) = f0(4, x, y_bot);
        f(5, x, y_bot) = f0(7, x, y_bot);
        f(6, x, y_bot) = f0(8, x, y_bot);
    }

    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            if (obstacle_mask(x, y)) {
                for (int i = 0; i < 9; ++i) {
                    f(i, x, y) = f0(opp[i], x, y);
                }
            }
        }
    }
}

void D2Q9_2D::macroscopic() {
    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            double rho_sum = 0.0, ux_sum = 0.0, uy_sum = 0.0;

            for (int i = 0; i < 9; ++i) {
                double fi = f(i, x, y);
                rho_sum += fi;
                ux_sum += fi * cx[i];
                uy_sum += fi * cy[i];
            }

            rho(x, y) = rho_sum;

            if (rho_sum != 0.0) {
                vx(x, y) = ux_sum / rho_sum + gravity * dt;
                vy(x, y) = uy_sum / rho_sum;
            } else {
                vx(x, y) = 0.0;
                vy(x, y) = 0.0;
            }
        }
    }
}

void D2Q9_2D::step() {
    macroscopic();
    equilibrium_f();
    collision();
    stream();
    bounce_boundaries();
}

void D2Q9_2D::run() {
    try {
        for (int t = 0; t < 1000; ++t) {
            step();
            write_field_to_file(rho, "output/rho", t);
            write_field_to_file(vx, "output/vx", t);
            write_field_to_file(vy, "output/vy", t);

            if (t % 100 == 0)
                std::cout << "Step: " << t << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Simulation error: " << e.what() << std::endl;
    }
}