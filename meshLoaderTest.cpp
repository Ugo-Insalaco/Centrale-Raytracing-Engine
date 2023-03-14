#include "meshLoader.cpp"
using namespace std;

int main(){
    TriangleMesh mesh("./Peachs Castle Exterior", "Peaches Castle");
    mesh.materials.at("Shape_131").texture.at(0).print();
    mesh.materials.at("Shape_131").texture.at(25).print();
    cout << mesh.indices.at(0).group << endl;
    return 0;
}