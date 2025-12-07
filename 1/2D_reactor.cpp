#include <2D_reactor.h>

reactor::reactor(int N_A, int N_B, double density_a, double density_b, double cutoff_radius): N_A(N_A), N_B(N_B), N(N_A + N_B), cutoff(cutoff_radius){
    cc=0;
    double total_density = (density_a+density_b);
    random_box box(N_A + N_B, total_density, 2, cutoff_radius);
    box_size=box.getBox();
    std::cout<<box_size<<std::endl;
    atoms=box.getItems();
    ca=N_A/(box_size*box_size);
    cb=N_B/(box_size*box_size);
}

void reactor::mark_a_b(){
    //g through box items and randomly mark them, keep indicies of what is A and what is B
    std::vector<int> indices(N_A + N_B); //it is assumed that all Na and Nb elements are placed
    for (int i = 0; i < N_A + N_B; ++i) {
        indices[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    for (int i = 0; i < N_A; ++i) {
        markers[indices[i]] = 'A'; 
    }
    for (int i = N_A; i < N_A + N_B; ++i) {
        markers[indices[i]] = 'B'; 
    }

}


void reactor::ideal_gas_collisions(){
    std::vector<bool> reacted(N, false);

    for (int idx=0; idx<N;idx++){
        if (reacted[idx]) continue;

        std::vector<double>& current_atom=atoms[idx];
        for (int jdx = 0; jdx < N; jdx++) {
            if (idx == jdx || reacted[jdx]) continue;

            std::vector<double>& atom = atoms[jdx];

            double dist=distance(current_atom,atom,box_size);

            
            if (dist<cutoff){
                bool current_in_a = markers.count(idx) && markers[idx] == 'A';
                bool other_in_b = markers.count(jdx) && markers[jdx] == 'B';

                if (current_in_a && other_in_b){
                    

                    std::vector<double> new_atom(2);  
                    new_atom[0] = (current_atom[0] + atom[0]) / 2.0;  
                    new_atom[1] = (current_atom[1] + atom[1]) / 2.0;
                    atoms[idx]=new_atom; //this is shit because i am keeping half empties in memory, alternative approach would be to have A and B and C seperate containters (might work better)
                    markers[idx] = 'C';
                    markers[jdx] = ' ';
                    reacted[idx] = true;  
                    reacted[jdx] = true;
                    break;
                }
            }            
        }
    }
}

void reactor::get_concentration(){
    int countA = 0, countB = 0,countC = 0;

    for ( auto& pair : markers) {
        if (pair.second == 'A') {
            countA++;
        } 
        else if (pair.second == 'B') {
            countB++;
        } 
        else if (pair.second == 'C') {
            countC++;
        }
    }
    ca=countA/(box_size*box_size);
    cb=countB/(box_size*box_size);
    cc=countC/(box_size*box_size);
}

void reactor::ideal_gas(){//atoms have no volume themselves
    //move all of them
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disr(0, space_step);
    std::uniform_real_distribution<> disfi(0, 2*M_PI);
    //for (auto& item : atoms){
    for (int idx = 0; idx < atoms.size(); idx++){
    if (markers[idx] != ' ') { 
        std::vector<double>& item = atoms[idx];
        double r=disr(gen);
        double fi=disfi(gen);
        item[0]+=r*std::cos(fi);
        item[1]+=r*std::sin(fi);

        item[0] -= box_size * std::floor((item[0] + box_size / 2.0) / box_size);
        item[1] -= box_size * std::floor((item[1] + box_size / 2.0) / box_size);
    }}
    ideal_gas_collisions();
}

void reactor::print_to_file(std::string path){
    std::ofstream outFile(path);

    for (size_t i = 0; i < atoms.size(); ++i) {
        char letter = markers[i]; 
        if (letter!=' '){
            const auto& pair = atoms[i]; 
            double x = pair[0];
            double y = pair[1];
            outFile << letter << "\t" << x << "\t" << y << std::endl;
        }
       
   }

   
   outFile.close();
}

void reactor::run(int steps, double space){
    space_step=space;
    mark_a_b();
    std::ofstream concfile("output/concnetration.txt");
    for (int i=0;i<=steps;i++){
        std::string filepath = "output/times/" + std::to_string(i) + ".txt";
        print_to_file(filepath);
        get_concentration();
        concfile << ca<<"\t"<<cb<<"\t"<<cc<<"\t"<< std::endl;
        if (i==steps){break;}
        ideal_gas();
    }
    concfile.close();
}

void reactor::run2(int steps, double space,std::string path){
    space_step=space;
    mark_a_b();
    std::ofstream concfile(path);
    for (int i=0;i<=steps;i++){
        std::string filepath = "output/times/" + std::to_string(i) + ".txt";
        get_concentration();
        concfile << ca<<"\t"<<cb<<"\t"<<cc<<"\t"<< std::endl;
        if (i==steps){break;}
        ideal_gas();
    }
    concfile.close();
}


/////////////////////////////////////////////////////////////////////////////////////////////////
porous_reactor::porous_reactor(int N_A, int N_B, int N_P, double density_a, double density_b,double cutoff_radius):N_A(N_A), N_B(N_B), N(N_A + N_B), cutoff(cutoff_radius){
    cc=0;
    double total_density = (density_a+density_b); //it needs to be so all have the same starting no matter how many P is added
    random_box box(N_A + N_B + N_P, total_density, 2, 1);//pores have radius 0.5!
    box_size=box.getBox();
    
    std::vector<std::vector<double>> atoms=box.getItems();
    ca=N_A/(box_size*box_size);
    cb=N_B/(box_size*box_size);
    
    mark_a_b(N_P,atoms);
}

void porous_reactor::mark_a_b(int N_P,std::vector<std::vector<double>> atoms){
    A.resize(N_A, std::vector<double>(2));
    B.resize(N_B, std::vector<double>(2));
    P.resize(N_P, std::vector<double>(2));
    std::vector<int> indices(N_A + N_B + N_P); 
    for (int i = 0; i < N_A + N_B + N_P; ++i) {
        indices[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    for (int i = 0; i < N_A; ++i) {
        A[i]=atoms[indices[i]]; 
    }

    for (int i = N_A; i < N_A + N_B; ++i) {
        B[i - N_A]=atoms[indices[i]];  
    }

    for (int i = N_A+N_B; i < N_A + N_B + N_P; ++i) {
        P[i - N_A - N_B]=atoms[indices[i]];  
    }
    std::cout<<A.size()<<std::endl;
}

bool porous_reactor::is_in_pore(std::vector<double>& partcile){
    for (auto& pore : P){
        double dist=distance(partcile,pore,box_size);
        if (dist<0.5){
            return true;
        }
    }
    return false;
}

void porous_reactor::ideal_gas(){//atoms have no volume themselves
    //move all of them
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disr(0, space_step);
    std::uniform_real_distribution<> disfi(0, 2*M_PI);
    //these 3 can run in parallel
    for (int ia = 0; ia < A.size(); ia++){
    
        std::vector<double>& item = A[ia];
        std::vector<double> new_item(2);
        double r=disr(gen);
        double fi=disfi(gen);
        new_item[0]=item[0]+r*std::cos(fi);
        new_item[1]=item[1]+r*std::sin(fi);

        new_item[1] -= box_size * std::floor((new_item[0] + box_size / 2.0) / box_size);
        new_item[1] -= box_size * std::floor((new_item[1] + box_size / 2.0) / box_size);
       
        if (!is_in_pore(new_item)){
            item[0]= new_item[0];
            item[1]= new_item[1];
        }
    }
    for (int ib = 0; ib < B.size(); ib++){
    
        std::vector<double>& item = B[ib];
        std::vector<double> new_item(2);
        double r=disr(gen);
        double fi=disfi(gen);
        new_item[0]=item[0]+r*std::cos(fi);
        new_item[1]=item[1]+r*std::sin(fi);

        new_item[1] -= box_size * std::floor((new_item[0] + box_size / 2.0) / box_size);
        new_item[1] -= box_size * std::floor((new_item[1] + box_size / 2.0) / box_size);
        
        if (!is_in_pore(new_item)){
            item[0]= new_item[0];
            item[1]= new_item[1];
        }
    }
    for (int ic = 0; ic < C.size(); ic++){
    
        std::vector<double>& item = C[ic];
        std::vector<double> new_item(2);
        double r=disr(gen);
        double fi=disfi(gen);
        new_item[0]=item[0]+r*std::cos(fi);
        new_item[1]=item[1]+r*std::sin(fi);

        new_item[0] -= box_size * std::floor((new_item[0] + box_size / 2.0) / box_size);
        new_item[1] -= box_size * std::floor((new_item[1] + box_size / 2.0) / box_size);
        
        if (!is_in_pore(new_item)){
            item[0]= new_item[0];
            item[1]= new_item[1];
        }
    }
    ideal_gas_collisions();
}

void porous_reactor::ideal_gas_collisions(){
    std::vector<bool> reactedA(N_A, false);
    std::vector<bool> reactedB(N_B, false);
    for (size_t ia=0;ia<A.size();++ia){
        if (reactedA[ia]==true) continue;
        std::vector<double>& a_current= A[ia];
        for (size_t ib=0;ib<B.size();++ib){
            if (reactedB[ib]==true) continue;
            std::vector<double>& b_current= B[ib];
            double dist=distance(a_current,b_current,box_size);
            if (dist<cutoff){
                std::vector<double> new_atom(2);  
                new_atom[0] = (a_current[0] + b_current[0]) / 2.0;  
                new_atom[1] = (a_current[1] + b_current[1]) / 2.0;
                if (!is_in_pore(new_atom)){
                    reactedA[ia]=true;
                    reactedB[ib]=true;
                    C.push_back(new_atom);
                    break;
                }

            }
        }
    }

    auto removeA = std::remove_if(A.begin(), A.end(), [&reactedA, index = 0](const std::vector<double>&) mutable {
        return reactedA[index++];
    });
    A.erase(removeA, A.end());

    auto removeB = std::remove_if(B.begin(), B.end(), [&reactedB, index = 0](const std::vector<double>&) mutable {
        return reactedB[index++];
    });
    B.erase(removeB, B.end());
}

void porous_reactor::run(int steps, double space){
    space_step=space;
    std::ofstream concfile("output/porous/concnetration.txt");
    for (int i=0;i<=steps;i++){
        std::string filepath = "output/porous/times/" + std::to_string(i) + ".txt";
        print_to_file(filepath);
        get_concentration();
        concfile << ca<<"\t"<<cb<<"\t"<<cc<<"\t"<< std::endl;
        if (i==steps){break;}
        
        ideal_gas();
    }
    concfile.close();
}

void porous_reactor::run2(int steps, double space,std::string path){
    space_step=space;
    std::ofstream concfile(path);
    for (int i=0;i<=steps;i++){
        
        concfile << ca<<"\t"<<cb<<"\t"<<cc<<"\t"<< std::endl;
        if (i==steps){break;}
        ideal_gas();
        get_concentration();
    }
    concfile.close();
}

void porous_reactor::get_concentration(){
   
    ca=A.size()/(box_size*box_size);
    cb=B.size()/(box_size*box_size);
    cc=C.size()/(box_size*box_size);
}



void porous_reactor::print_to_file(std::string directory){

    std::ofstream outFile(directory);
    for (size_t i = 0; i < A.size(); ++i) {
        const auto& pair = A[i]; 
        outFile <<"A"<<"\t" << pair[0] << "\t" << pair[1] << std::endl;
        }
    for (size_t i = 0; i < B.size(); ++i) {
        const auto& pair = B[i]; 
        outFile <<"B"<<"\t" << pair[0] << "\t" << pair[1] << std::endl;
    }
    for (size_t i = 0; i < C.size(); ++i) {
        const auto& pair = C[i]; 
        outFile <<"C"<<"\t" << pair[0] << "\t" << pair[1] << std::endl;
    }
    for (size_t i = 0; i < P.size(); ++i) {
        const auto& pair = P[i]; 
        outFile <<"P"<<"\t" << pair[0] << "\t" << pair[1] << std::endl;
    }
    outFile.close();
}