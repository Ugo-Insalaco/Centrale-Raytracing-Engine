#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

#include "triangleIndices.h"
#include "material.h"
#include <map> 

using namespace std;

class TriangleMesh {
public:
    ~TriangleMesh() {}
    TriangleMesh(string shapeName);
    void readMTL(string mtl);
    void readOBJ(string obj);
    void scale(double scale_ratio);
    void translate(Vector translation);
    Vector relativeTranslate(Vector translation);
    void rotate(Vector center, Vector rotation);
    void relativeRotate(Vector center, Vector rotation);
    void loadTextures(string path);

    string path;
    string fileName;
    vector<TriangleIndices> indices;
    map<string, Material*> materials;
    vector<Vector> vertices;
    vector<Vector> normals;
    vector<Vector> uvs;
    vector<Vector> vertexcolors;
    Vector rotation;
};      
#endif