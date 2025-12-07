#include <heat_conduction.h>

rectangular_rod::rectangular_rod(int N, double x0, double xn, double y0, double yn) : N(N), x0(x0), xn(xn), y0(y0), yn(yn), temperatures(N,N), starting_temeperatures(N,N){

}

void rectangular_rod::boundary_values(double h){
    for (int i=0;i<N;i++){
        starting_temeperatures(i,0)-=x0;
    }
    for (int i=0;i<N;i++){
        starting_temeperatures(i,N-1)-=xn;
    }
    for (int i=0;i<N;i++){
        starting_temeperatures(0,i)-=y0;
    }
    for (int i=0;i<N;i++){
        starting_temeperatures(N-1,i)-=yn;
    }
}

void rectangular_rod::boundary_values_neuman(double h){
    for (int i=0;i<N;i++){
        starting_temeperatures(i,0)-=x0;
    }
    for (int i=0;i<N;i++){
        starting_temeperatures(i,N-1)-=xn;
    } 
    for (int i=0;i<N;i++){
        starting_temeperatures(0,i)-=y0;
    }
    for (int i=0;i<N;i++){
        starting_temeperatures(N-1,i)-=yn;
    }
}



void rectangular_rod::solver_laplace(double h){
    //Make A
    csr A;
    A.tridiagonal(N,1,-4,1);
    //Make I
    csr I;
    I.tridiagonal(N,0,1,0);

    csr D;
    D.block_tridiagonal(A,I);
    //dealloacate A, I

    //D.print_matrix(N*N,N*N); this is ok

    boundary_values(h);
    //std::cout<<starting_temeperatures;
    

    std::vector<double> temps =GC_solver(D,starting_temeperatures.get_data());
    temperatures.set(temps);
    //std::cout<<temperatures<<std::endl; //this is not right, problem with either solver or the way data is habdled
    //dealoacate temps
}

void rectangular_rod::solver_poisson(double h, double source){
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            starting_temeperatures(i,j)=-source*h*h;
        }
    }


    //Make A
    csr A;
    A.tridiagonal(N,1,-4,1);
    //Make I
    csr I;
    I.tridiagonal(N,0,1,0);

    csr D;
    D.block_tridiagonal(A,I);
    //dealloacate A, I

    //D.print_matrix(N*N,N*N); this is ok

    boundary_values(h);
    //std::cout<<starting_temeperatures;
    

    std::vector<double> temps =GC_solver(D,starting_temeperatures.get_data());
    temperatures.set(temps);
    //std::cout<<temperatures<<std::endl; //this is not right, problem with either solver or the way data is habdled
    //dealoacate temps
}

void rectangular_rod::solver_laplace_neuman(double h){
    //Make A
    csr A;
    A.tridiagonal(N,1,-4,1);
    A.modify_last_diagonal(-3);
    //Make I
    csr I;
    I.tridiagonal(N,0,1,0);

    csr D;
    D.block_tridiagonal(A,I);
    //dealloacate A, I

    /*
    for (int i = 0; i < N; ++i) {
        int row = i * N + (N - 1); // Row index for (i, N-1) in row-major order
        // Find and modify the diagonal entry
        for (int j = D.get_row_at(row); j < D.get_row_at(row + 1); ++j) {
            if (D.get_col_at(j) == row) { // Diagonal entry
                D.set_value_at(j, -3);    // Modify value directly
                break;
            }
        }
    }
    */

    //D.print_matrix(N*N,N*N); //this is ok

    boundary_values_neuman(h);
    //std::cout<<starting_temeperatures;
    

    std::vector<double> temps =GC_solver(D,starting_temeperatures.get_data());
    temperatures.set(temps);
    //std::cout<<temperatures<<std::endl; //this is not right, problem with either solver or the way data is habdled
    //dealoacate temps
}

void rectangular_rod::solver_poisson_neuman(double h, double source){
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            starting_temeperatures(i,j)=-source*h*h;
        }
    }
    //Make A
    csr A;
    A.tridiagonal(N,1,-4,1);
    A.modify_last_diagonal(-3);
    //Make I
    csr I;
    I.tridiagonal(N,0,1,0);

    csr D;
    D.block_tridiagonal(A,I);
    //dealloacate A, I

    /*
    for (int i = 0; i < N; ++i) {
        int row = i * N + (N - 1); // Row index for (i, N-1) in row-major order
        // Find and modify the diagonal entry
        for (int j = D.get_row_at(row); j < D.get_row_at(row + 1); ++j) {
            if (D.get_col_at(j) == row) { // Diagonal entry
                D.set_value_at(j, -3);    // Modify value directly
                break;
            }
        }
    }
    */

    //D.print_matrix(N*N,N*N); //this is ok

    boundary_values_neuman(h);
    //std::cout<<starting_temeperatures;
    

    std::vector<double> temps =GC_solver(D,starting_temeperatures.get_data());
    temperatures.set(temps);
    //std::cout<<temperatures<<std::endl; //this is not right, problem with either solver or the way data is habdled
    //dealoacate temps
}

void rectangular_rod::print_to_file(std::string path){
    std::ofstream outFile(path);
    outFile << temperatures;
    outFile.close();
}

