#include <FFT.h>
#include <iostream>
#include <fstream>

int main(){

    std::ifstream file("C:/Users/domin/OneDrive/Namizje/MKS/2/output/duffing/6.txt"); 
    if (!file.is_open()) {
        std::cerr << "Failed to open file.\n";
        return 1;
    }

    std::vector<std::complex<double>> x_values;
    double t, x, dx;

    while (file >> t >> x >> dx) {
        x_values.push_back(x);
    }

    file.close();
    dfft(x_values);
    std::ofstream outfile("6.txt");
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return 1;
    }
    for (const auto& x : x_values) {
        outfile << real(x) << "\t" << imag(x) << std::endl;
    }

    outfile.close();
}