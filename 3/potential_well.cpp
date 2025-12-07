#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#define pi 3.141592653589

double  dpsi(double x, double psi, double z){
    return z;
}

double dz(double x, double psi, double z, double E){
    return -2.0*E*psi;
}

double krajisce(double E, bool odd, double& normal){
    double x=0;
    double psi;
    double z;

    if (odd==0){
        psi=0;
        z=1;
    }
    else{
        psi=1;
        z=0;
    }

    double h=0.001;
    std::string path="jama.txt";
    std::ofstream outFile(path);
    

    double integral=psi*psi*h*0.5;

    while (x<=1){
        outFile << x << "\t" << psi<<"\t"<<z << std::endl;

        if (odd==0){
            outFile << -x << "\t" << -psi<<"\t"<<z << std::endl;
        }
        else{
            outFile << -x << "\t" << psi<<"\t"<<z << std::endl;
        }

        double k1=h*dpsi(x,psi,z);
        double l1=h*dz(x,psi,z,E);
        double k2=h*dpsi(x+h/2.0,psi+k1/2.0,z+l1/2.0);
        double l2=h*dz(x+h/2.0,psi+k1/2.0,z+l1/2.0,E);
        double k3=h*dpsi(x+h/2.0,psi+k2/2.0,z+l2/2.0);
        double l3=h*dz(x+h/2.0,psi+k2/2.0,z+l2/2.0,E);
        double k4=h*dpsi(x+h,psi+k3,z+l3);
        double l4=h*dz(x+h,psi+k3,z+l3,E);

        x+=h;
        psi+=1/6.0 * (k1+2*k2+2*k3+k4);
        z+=1/6.0 * (l1+2*l2+2*l3+l4);

        integral+=psi*psi*h;

    }
    integral-=psi*psi*h*0.5;
    integral*=2;
    outFile.close();

    normal=integral;

    return psi;

}



int main(){
    int odd;
    std::cout<<"Odd or even [0/1]"<<std::endl;

    std::cin>>odd;

    double El,Er;
    int n;
    std::cout<<"energy state"<<std::endl;
    std::cin>> n;

    El=pi*pi*n*n/8.0 - 0.1;
    Er=pi*pi*n*n/8.0 + 0.1;

    double tol=1e-5;
    double c;
    double normal=1;
    while (std::abs(El-Er)>tol){
        c=(El+Er)/2.0;
        if (krajisce(El,odd,normal)*krajisce(c,odd,normal)<0){
            Er=c;
        }
        else{
            El=c;
        }
    }
    std::cout<<c<<std::endl;
    std::cout<<normal<<std::endl;
}

