using namespace std;
#include <random>
#include <iostream>
#include <cmath>

const double pi = 3.141592654;
double Gaussian(double x){
    return 1/sqrt(2*pi)*exp(-pow(x, 2)/2);
}

double Integrande(double x, double y, double z){
    if(x <= -pi/2 || x >= pi/2) x = 0;
    if(y <= -pi/2 || y >= pi/2) y = 0;
    if(z <= -pi/2 || z >= pi/2) z = 0;
    return pow(cos(x+y+z), 5);
}

int main(){
    int N=500;
    double s = 0;
    for (int i = 0; i < N; i++){
        double x1, y1, z1, x2, y2, z2;
        double r1, r2;

        r1 = (float) rand()/RAND_MAX, r2 = (float) rand()/RAND_MAX;
        x1 = sqrt(-2*log(r1))*cos(2*pi*r2), x2 = sqrt(-2*log(r1))*sin(2*pi*r2);
        // cout << x1 << " " << x2 << endl;
        r1 = (float) rand()/RAND_MAX, r2 = (float) rand()/RAND_MAX;
        y1 = sqrt(-2*log(r1))*cos(2*pi*r2), y2 = sqrt(-2*log(r1))*sin(2*pi*r2);

        r1 = (float) rand()/RAND_MAX, r2 = (float) rand()/RAND_MAX;
        z1 = sqrt(-2*log(r1))*cos(2*pi*r2), z2 = sqrt(-2*log(r1))*sin(2*pi*r2);

        s += Integrande(x1, y1, z1)/(Gaussian(x1)*Gaussian(y1)*Gaussian(z1));
        cout << Integrande(x1, y1, z1)/(Gaussian(x1)*Gaussian(y1)*Gaussian(z1)) << endl;
        s += Integrande(x2, y2, z2)/(Gaussian(x2)*Gaussian(y2)*Gaussian(z2));
    }
    cout << s/2/N;
    int t;
    cin >> t;
    return 0;
}