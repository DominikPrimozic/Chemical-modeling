#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <chrono>


std::complex<double> g(std::complex<double> z){
    return std::pow(z,3)-std::complex<double>(1.0, 0.0);
}

std::complex<double> dg(std::complex<double> z){
    return 3.0*std::pow(z,2);
}

std::pair<std::complex<double>, int> newthon(std::complex<double> z0, double tol=1e-9, int maxit=1000){
    std::complex<double> z;
    int iterations = 0;
    for (int it=0;it<maxit;it++){
        z= z0-g(z0)/dg(z0);
        if  (std::abs(z-z0)<tol){break;}
        z0=z;
        iterations++;
    }
    return make_pair(z, iterations);
}


int rootindex(std::complex<double> z, std::vector<std::complex<double>>& roots){
    for (int i=0;i<roots.size();i++){
        if (abs(roots[i]-z)<1e-5){return i;}
    }
    #pragma omp critical
    {
        roots.push_back(z);
    }
    return roots.size()-1;
}

void runer(int res, double dx, double dy, double xm, double ym, int num_threads){
    omp_set_num_threads(num_threads);
    std::vector<std::ostringstream> thread_buffers(3);

    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<std::complex<double>> roots;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<std::ostringstream> local_buffers(3);

        #pragma omp for
        for (int k = 0; k < res * res; ++k) {
            int i = k / res;    // Compute the row index
            int j = k % res;    // Compute the column index
            double x = xm + i * dx;
            double y = ym + j * dy;

            std::complex<double> z0(x, y);
            std::pair<std::complex<double>, int> z = newthon(z0);
            int root = rootindex(z.first, roots); 
            int iterations = z.second;

            // Accumulate results into thread-local buffers
            local_buffers[root] << x << "\t" << y << "\t" << iterations << "\n";
        }

        #pragma omp critical
        {
            for (int i = 0; i < 3; ++i) {
                thread_buffers[i] << local_buffers[i].str();
            }
        }
    }

    
    for (int i = 0; i < 3; ++i) {
        std::ofstream file("root" + std::to_string(i) + ".txt");
        file << thread_buffers[i].str();
        file.close();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;

    std::cout << "Total time taken: " << total_time.count() << " seconds" << std::endl;
}
void runer2(int res, double dx, double dy, double xm, double ym, int num_threads){
    omp_set_num_threads(num_threads);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<std::complex<double>> roots;

    #pragma omp parallel //or just go #pragma omp parallel for //on the loop
    {
        #pragma omp for
        for (int k = 0; k < res * res; ++k) {
            int i = k / res;    // Compute the row index
            int j = k % res;    // Compute the column index
            double x = xm + i * dx;
            double y = ym + j * dy;

            std::complex<double> z0(x, y);
            std::pair<std::complex<double>, int> z = newthon(z0);
            int root = rootindex(z.first, roots); //for real speed comparison i should have a vector with roots and just compare them, cause this can easily go too many

        }
    }

    

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;

    std::cout << "Total time taken: " << total_time.count() << " seconds" << std::endl;

    std::cout << "Identified roots: " << std::endl;
    for (int i = 0; i < roots.size(); ++i) {
        std::cout << "Root " << i << ": " << roots[i] << std::endl;
    }
}

int main(){
    std::cout<<omp_get_max_threads()<<std::endl;
    int num_threads = 4;  // Adjust this based on the number of processors you want to use
    omp_set_num_threads(num_threads);

    std::vector<std::complex<double>> roots;

    double xm=-2,xp=2,ym=-2,yp=2;
    int res=800;
    double dx=(xp-xm)/(res-1);
    double dy=(yp-ym)/(res-1);
    /*
    std::ofstream file_root0("root0.txt");  // Will dynamically use these files based on root discovery
    std::ofstream file_root1("root1.txt");
    std::ofstream file_root2("root2.txt");
    auto start_time = std::chrono::high_resolution_clock::now();

    
    
    
    #pragma omp for 
    for (int k = 0; k <res * res; ++k) {
        int i = k / res;    // Compute the row index
        int j = k % res;    // Compute the column index
        double x = xm + i * dx;
        double y = ym + j * dy;

        std::complex<double> z0(x,y);

        std::pair<std::complex<double>, int> z=newthon(z0);

        // int root = classify_or_add_root(z.first, roots);
        int root= rootindex(z.first, roots);
        int iterations = z.second; 

            //write to file coordinates from which it goes to the zero
        #pragma omp critical  // Ensure one thread writes at a time
        {
            if (root == 0) {
                file_root0 << x << "\t" << y <<"\t\t"<<iterations<< "\n"; 
            } else if (root == 1) {
                file_root1 << x << "\t" << y <<"\t\t"<<iterations<< "\n";  
            } else if (root == 2) {
                file_root2 << x << "\t" << y <<"\t\t"<<iterations<< "\n";  
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;

    
    file_root0.close();
    file_root1.close();
    file_root2.close();

    std::cout << "Identified roots: " << std::endl;
    for (int i = 0; i < roots.size(); ++i) {
        std::cout << "Root " << i << ": " << roots[i] << std::endl;
    }

    std::cout << "Total time taken: " << total_time.count() << " seconds" << std::endl;*/

    runer2(res,dx,dy,xm,ym,1); 
    // runer(res,dx,dy,xm,ym,2); 
     // runer(res,dx,dy,xm,ym,3); 
     //  runer(res,dx,dy,xm,ym,4); 
      //  runer(res,dx,dy,xm,ym,5); 
         runer2(res,dx,dy,xm,ym,6); 
          runer2(res,dx,dy,xm,ym,7); 
           runer(res,dx,dy,xm,ym,8); 
            runer2(res,dx,dy,xm,ym,9); 
            // runer(res,dx,dy,xm,ym,10); 
             // runer(res,dx,dy,xm,ym,11); 
               //runer(res,dx,dy,xm,ym,12); //it is not faster :(, i think the slow part is file writing 
    
}