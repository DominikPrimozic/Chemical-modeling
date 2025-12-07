#include <random_walk.h>
#include <regex>

double distance_nonPBS(const std::vector<double>& coord1, const std::vector<double>& coord2) {
    double dist_sq = 0.0;
    for (size_t i = 0; i < coord1.size(); ++i) {
        double delta = coord1[i] - coord2[i];
        dist_sq += delta * delta;
    }
    return sqrt(dist_sq);
}

bool walker::isPositionValid(const std::vector<double>& pos) {
    for (const auto& c : obstacles) {
        if (distance_nonPBS(pos, c) <= ball_radius) { 
            return false;
        }
    }
    return true;
}

std::vector<double> walker::adjustStartPosition() {
    std::vector<double> startPos = {0, 0};
    
    if (!isPositionValid(startPos)) {
        double shift = 1e-3; 
        int max_attempts = 5000;

        for (int i = 0; i < max_attempts; ++i) {
            startPos[0] += shift;
            startPos[1] += shift;  

            if (isPositionValid(startPos)) {
                return startPos;
            }
        }
        
    }

    return startPos;
}


walker::walker(int N, double density, double limit_size, int reproducible) : ball_radius(limit_size/2.0) {

    path_path="output/random_walk/path/" +  std::to_string(N) + "_" + std::to_string(density)+ "_" + std::to_string(reproducible) + ".txt";
    distance_path="output/random_walk/distance/" + std::to_string(N) + "_" + std::to_string(density) + "_" + std::to_string(reproducible)+ ".txt";
    steps_path="output/random_walk/steps/" + std::to_string(N) + "_" + std::to_string(density) + "_" + std::to_string(reproducible)+ ".txt";
    configuration_path="output/random_walk/configuration/" +  std::to_string(N) + "_" + std::to_string(density) +"_" + std::to_string(reproducible)+ ".txt";

    random_box placement(N,density,2, limit_size);
    obstacles=placement.getItems();
    box_size=placement.getBox();
    std::cout<<box_size<<std::endl;
    std::cout<<ball_radius<<std::endl;
    obstacles_print();
    //would prefer to dealoacet placement here
    //position={0,0};
    position = adjustStartPosition();
    origin=position;
    statistial_position=position;
}

void walker::step(double step_size){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 360);
    bool ismove=false;

    double x, y;

    while (!ismove){ 
        bool couldnotmove=false;
        double move=dis(gen);
               x=step_size*std::cos(move*pi/180.0);
               y=step_size*std::sin(move*pi/180.0);

        for (const auto& c : obstacles){
            if (distance_nonPBS({position[0] + x,position[1] + y},c) <=ball_radius ){ //limit_size/2.0
                couldnotmove=true;
                failed++;
                break;
            }
        }
        if (!couldnotmove){
            if (position[0] + x >= -box_size/2.0 && position[0] + x <= box_size/2.0 && position[1] + y >= -box_size/2.0 && position[1] + y <= box_size/2.0){
                ismove=true;
            }
            else failed++;
        }
        if (failed>10000) break; //should error here though, std::cerr << "Failed too many times";
            
    }
    position[0]+=x;
    position[1]+=y;
}

void walker::walk(int steps,double step_size){
    int current_step=0;
    distance_from_0=0;
    std::ofstream path_output(path_path);
    std::ofstream distance_output(distance_path);
    std::ofstream steps_output(steps_path);
    path_output << position[0]<<"\t"<< position[1]<<std::endl;
    distance_output<<distance_from_0<<std::endl;
    steps_output<<steps<<"\n"<<std::endl;
    while (current_step<steps){
        step(step_size);
        distance_from_origin();
        path_output << position[0]<<"\t"<< position[1]<<std::endl;
        distance_output<<distance_from_0<<std::endl;
        steps_output<<failed<<std::endl;
        current_step++;
        failed=0;
    }
    path_output.close();
    distance_output.close();
    steps_output.close();
}

void walker::obstacles_print(){
    std::ofstream outfile(configuration_path);

    for (int c=0;c<obstacles.size();c++) {
        for (int d=0;d<obstacles[c].size();d++) {
            outfile << obstacles[c][d] << "\t";
        }
        outfile << std::endl;
    }
    outfile.close();
}

void walker::distance_from_origin(){
    distance_from_0=0;
    for (int i=0;i<position.size();i++){
        distance_from_0+=(origin[i]-position[i])*(origin[i]-position[i]);
    }
    distance_from_0=std::sqrt(distance_from_0);
}

std::string updateFilename(const std::string& original, int rep) { //this i took from chatgpt
    std::regex num_pattern(R"(_\d+\.txt$)"); // Matches _123.txt at the end
    std::string new_suffix = "_" + std::to_string(rep) + ".txt";

    if (std::regex_search(original, num_pattern)) {
        // If a number is already appended, replace it
        return std::regex_replace(original, num_pattern, new_suffix);
    } else {
        // Otherwise, insert the number before ".txt"
        size_t pos = original.rfind(".txt");
        if (pos != std::string::npos) {
            return original.substr(0, pos) + new_suffix;
        }
    }
    return original; // Fallback if ".txt" is not found
}

void walker::multiple_walks(int walkN,int steps,double step_size){
    for (int rep=0;rep<walkN;rep++){
        position = statistial_position;
        path_path=updateFilename(path_path, rep); 
        steps_path=updateFilename(steps_path, rep); 
        configuration_path=updateFilename(configuration_path, rep); 
        distance_path=updateFilename(distance_path, rep); 
        walk(steps,step_size);
        
    }
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string updateFilename2(const std::string& original, double rep) {
    std::regex num_pattern(R"(_\d+\.txt$)"); // Matches _123.txt at the end
    std::string new_suffix = "_" + std::to_string(rep) + "_"+ "0" +".txt";

    if (std::regex_search(original, num_pattern)) {
        // If a number is already appended, replace it
        return std::regex_replace(original, num_pattern, new_suffix);
    } else {
        // Otherwise, insert the number before ".txt"
        size_t pos = original.rfind(".txt");
        if (pos != std::string::npos) {
            return original.substr(0, pos) + new_suffix;
        }
    }
    return original; // Fallback if ".txt" is not found
}


std::string updateFilename3(const std::string& original, double new_value) {
    std::regex pattern(R"((_\d+\.\d+)_([\d\.]+)_\d+\.txt$)"); // Matches _X.Y_Z.W.txt
    std::smatch match;

    if (std::regex_search(original, match, pattern)) {
        std::string new_suffix = match[1].str() + "_" + std::to_string(new_value) + "_0.txt";
        return std::regex_replace(original, pattern, new_suffix);
    }

    return original; // Return unchanged if no match found
}




void walker::update_paths(double step_size){
    path_path=updateFilename2(path_path, step_size); 
    steps_path=updateFilename2(steps_path, step_size); 
    configuration_path=updateFilename2(configuration_path, step_size); 
    distance_path=updateFilename2(distance_path, step_size);
}

void walker::multiple_walks2(int walkN,int steps,double step_size){
    path_path=updateFilename3(path_path, step_size); 
    steps_path=updateFilename3(steps_path, step_size); 
    configuration_path=updateFilename3(configuration_path, step_size); 
    distance_path=updateFilename3(distance_path, step_size); 
    for (int rep=0;rep<walkN;rep++){
        position = statistial_position;
        path_path=updateFilename(path_path, rep); 
        steps_path=updateFilename(steps_path, rep); 
        configuration_path=updateFilename(configuration_path, rep); 
        distance_path=updateFilename(distance_path, rep); 
        walk(steps,step_size);
    
    }
}