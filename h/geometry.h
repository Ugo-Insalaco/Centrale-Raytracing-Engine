#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "enums.h"
class Geometry{
    public:
        virtual Vector intersection(Ray ray, float& t, Vector& normal) = 0;
        Vector position;
        Vector albedo;
        Texture texture;
        float refractionIndex;
};
#endif