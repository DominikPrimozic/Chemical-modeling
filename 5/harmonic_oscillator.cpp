#include <harmonic_oscillator.h>

double hermite(double x, int n){
    if (n==0) return 1;
    if (n==1) return 2*x;
    return 2*x*hermite(x,n-1) - 2*(n - 1)*hermite(x,n-2);
}

double factorial(int n){
    if (n==0) return 1;
    return factorial(n-1);
}

double x_n(double x, int n){
    double first= 1/std::sqrt(std::pow(2,n)*factorial(n)*std::sqrt(M_PI));
    double gaussian = std::exp(-x*x/2);
    double herm = hermite(x,n);
    return first*gaussian*herm;
}

quantum_harmonic::quantum_harmonic(int num): ap(num,num),am(num,num), H0(num,num){
    H0(0,0)=0.5;
    for (int n=1;n<num;n++){
        ap(n-1,n)=std::sqrt(n);
        am(n,n-1)=std::sqrt(n);
        H0(n,n)=0.5;
    }
}

matrix quantum_harmonic::perturbation_3(double alpha){
    return (alpha/(2*std::sqrt(2)))*(ap*ap*ap + 3*ap*ap*am + 3*ap*am*am + am*am*am);
}

void quantum_harmonic::solve_states(double alpha){
    matrix H=(am*ap) + H0 + perturbation_3(alpha);
    jacobi_solver js(H);
    js.solve();
    std::ofstream outfileV("eigenvectors.txt");
    outfileV<<js.get_eigenvectors()<<std::endl;
    outfileV.close();
    std::ofstream outfilel("eigenvalues.txt");
    outfilel<<js.get_eigenvalues()<<std::endl;
    outfilel.close();
}

void quantum_harmonic::get_states(double alpha, double L, int N){
    matrix H=(am*ap) + H0 + perturbation_3(alpha);
    jacobi_solver js(H);
    js.solve();
    std::ofstream outfileV("eigenvectors.txt");
    outfileV<<js.get_eigenvectors()<<std::endl;
    outfileV.close();
    std::ofstream outfilel("eigenvalues.txt");
    outfilel<<js.get_eigenvalues()<<std::endl;
    outfilel.close();

    double dx=2*L/N;
    const matrix& V = js.get_eigenvectors();
    for (int n=0;n<am.rows();n++){
        std::string path = "eigenvectors/" + std::to_string(n) + ".txt";
        std::ofstream outfile(path);
        double x0=-L;
        while (x0<=L){
            double psi=0;
            for (int nn=0;nn<V.rows();nn++){
                psi+=V(nn,n)*x_n(x0,nn);
            }
            outfile << x0 << "\t" << psi << std::endl;
            x0+=dx;
        }
        outfile.close();
    }
}
