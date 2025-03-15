#include <random_placer.h>



double distance(const std::vector<double>& coord1, const std::vector<double>& coord2, double l) {
    double dist_sq = 0.0;
    for (size_t i = 0; i < coord1.size(); ++i) {
        //double delta = fabs(coord1[i] - coord2[i]);
        double delta = coord1[i] - coord2[i];
       //delta = std::min(delta, l - delta);  
        delta = delta - std::floor((delta + l/2.0) / l) * l;
        dist_sq += delta * delta;
    }
    return sqrt(dist_sq);
}

random_box::random_box(int N, double density, int dimensionality,double limit_size) : dim(dimensionality), N(N){
    box_size=std::pow((N/density),1.0/dimensionality);
    place(limit_size);
    packing(density);
    //file_output();
}

void random_box::place(double limit_size){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-box_size/2.0, box_size/2.0);

    std::vector<double> coordinate(dim);
    for (int d = 0; d < dim; ++d) {
        coordinate[d] = dis(gen);
    }
    items.push_back(coordinate);

    int max_attempts=250;
    int attempts = 0;

    while (items.size()<N){
        if (attempts >= max_attempts) {
            items.clear();
            break;
        }
        for (int d = 0; d < dim; ++d) {
            coordinate[d] = dis(gen);
        }

        bool is_valid=true;
        for (std::vector<double> c : items){
            if (distance(coordinate, c, box_size)<limit_size){
                is_valid=false;
                break;
            }
        }
        if (is_valid){
            items.push_back(coordinate);
            attempts=0;
        }
        else{attempts++;}
    }
}

void random_box::packing(double density){
    double V = (N / density);
    double oneV;
    switch (dim) {
        case 1:
            oneV = 1;
            break;
        case 2:
            oneV = 3.141592653589 * 0.5 * 0.5;
            break;
        default:
            oneV = 4.0 / 3.0 * 3.141592653589 * 0.5 * 0.5 * 0.5;
            break;
    }

    packed=oneV*N/V;
}

void random_box::file_output(){
    std::ofstream outfile("output.txt");

    for (int c=0;c<items.size();c++) {
        for (int d=0;d<items[c].size();d++) {
            outfile << items[c][d] << "\t";
        }
        outfile << std::endl;
    }

    outfile <<"\n\n" <<std::endl;
    outfile << "Packing: " << packed << std::endl;
    outfile.close();
}

void random_box::file_output2(std::string path){
    std::ofstream outfile(path);

    for (int c=0;c<items.size();c++) {
        for (int d=0;d<items[c].size();d++) {
            outfile << items[c][d] << "\t";
        }
        outfile << std::endl;
    }

    outfile <<"\n\n" <<std::endl;
    outfile << "Packing: " << packed << std::endl;
    outfile.close();
}
bool random_box::file_output3(){
    return items.empty();
}
/*
random_box::random_box(int N, double density, int dimensionality) : dim(dimensionality), N(N) {
    box_size = std::pow((N / density), 1.0 / dimensionality);
    // Pre-allocate memory for N coordinates (no initialization to zero)
    items.resize(N); 
    place(); 
    file_output();  
}

void random_box::place() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-box_size / 2.0, box_size / 2.0);

    std::vector<double> coordinate(dim);  
    int max_attempts = 250;
    int attempts = 0;

    int placed_items = 0;  

    while (placed_items < N) {
        if (attempts >= max_attempts) {
            break;
        }

        
        for (int d = 0; d < dim; ++d) {
            coordinate[d] = dis(gen);
        }

        bool is_valid = true;

       
        for (int i = 0; i < placed_items; ++i) {
            if (distance(coordinate, items[i], box_size) < 1.0) {
                is_valid = false;
                break;
            }
        }

        
        if (is_valid) {
            items[placed_items] = coordinate;  
            placed_items++;  
            attempts = 0;  
        } else {
            attempts++; 
        }
    }
}
*/