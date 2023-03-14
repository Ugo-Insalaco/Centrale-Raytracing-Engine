#ifndef MATERIAL_H
#define MATERIAL_H
#include "vector.h"
#include <vector>
using namespace std;

class Material {
    public:
        Material();
        int width=0, height=0;
        Vector Kd;
        string KdMap="";
        vector<Vector> texture;
};
#endif