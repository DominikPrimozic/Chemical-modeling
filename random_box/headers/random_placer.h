#ifndef random_placement
#define random_placement

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random_placer.h>
#include <string>

class random_box{
    private:
    std::vector<std::vector<double>> items;
    double box_size;
    int dim;
    int N;
    double packed;

    void place(double limit_size);
    void packing(double density);

    public:
        random_box(int N, double density, int dimensionality, double limit_size=1);

        void file_output();
        void file_output2(std::string path);
        bool file_output3();

        std::vector<std::vector<double>> getItems() { return items; } 
        double getBox(){return box_size;}


};

double distance(const std::vector<double>& coord1, const std::vector<double>& coord2, double l);

#endif // random_placement