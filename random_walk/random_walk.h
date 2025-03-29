#ifndef twoD_walk
#define twoD_walk

#include <random_placer.h>
#include <vector>
#include <cmath>

#define pi 3.141592653589

class walker{
    private:
        std::vector<std::vector<double>> obstacles;
        std::vector<double> position, origin, statistial_position;
        int failed;
        double box_size;
        std::string path_path, distance_path, steps_path,configuration_path;
        double distance_from_0;
        double ball_radius;
        

        bool isPositionValid(const std::vector<double>& pos);
        std::vector<double> adjustStartPosition();
        

        void step(double step_size);
        void obstacles_print();
        void distance_from_origin();



    public:
        void walk(int steps,double step_size);
        walker(int N, double density, double limit_size=1,int reproducible=0);
        void multiple_walks(int walkN,int steps,double step_size);
        void multiple_walks2(int walkN,int steps,double step_size);

        void update_paths(double step_size);
        


};

#endif // twoD_walk
