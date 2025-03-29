#ifndef two_D_reactor
#define two_D_reactor

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random_placer.h>
#include <string>
#include <unordered_map>



#define M_PI           3.14159265358979323846


class reactor{ //ideal gas reactor
    private:
        int dimensionality=2;
        int N_A,N_B,N;
        double ca,cb,cc;
        double cutoff;
        double box_size;
        double space_step;
        std::vector<std::vector<double>> atoms;
        //std::vector<int> marked_a,marked_b,marked_c; //this is shit for memory and speed
        std::unordered_map<int, char> markers;


    public:
        reactor(int N_A, int N_B, double density_a, double density_b,double cutoff_radius);

        void run(int steps, double space);
        void run2(int steps, double space,std::string path);

        void mark_a_b();
        void get_concentration();

        void ideal_gas();
        void ideal_gas_collisions();
        
        void print_to_file(std::string path);



};

class porous_reactor{ //ideal gas reactor
    private:
        int dimensionality=2;
        int N_A,N_B,N;
        double ca,cb,cc;
        double cutoff;
        double box_size;
        double space_step;
        std::vector<std::vector<double>> A,B,C,P;


    public:
        porous_reactor(int N_A, int N_B, int N_P, double density_a, double density_b,double cutoff_radius);

        void run(int steps, double space);
        void run2(int steps, double space,std::string path);

        void mark_a_b(int N_P, std::vector<std::vector<double>> atoms);
        bool is_in_pore(std::vector<double>& partcile);
        void get_concentration();

        void ideal_gas();
        void ideal_gas_collisions();
        
        void print_to_file(std::string directory);



};



#endif // two_D_reactor

//Na/V, Nb/V