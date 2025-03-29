#ifndef probability_distributions
#define probability_distributions

#include <random>
#include <vector>
#include <omp.h>
#include <string>
#include <fstream>
#include <iostream>

# define PI           3.14159265358979323846

class maxwell_distribution{
    private:
        double k,m,T;
        double sigma;
        std::vector<int> x_bins;
        std::vector<int> v_bins;
        std::vector<int> E_bins;

    public:
    maxwell_distribution(int N, std::string path, double T, double m, double min_val,double max_val);
    void generate_distribution(int N, double min_val,double max_val,std::string path);
    void print_to_file(std::string path,int bins, std::vector<int>& bin_vector);

};

struct boxmuller_double {
    double z1;
    double z2;

    boxmuller_double(double u1, double u2);
}; //what is better: tuple, pair, structure?

int value_to_bin(double val, int bins, double max_val,double min_val);


#endif // probability_distributions
