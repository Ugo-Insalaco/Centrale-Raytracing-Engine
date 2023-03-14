#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H
#include "vector.h"
#include "ray.h"
class BoundingBox{
    public :
        BoundingBox();
        BoundingBox(Vector& bmin, Vector& bmax);
        bool intersection(Ray ray);
        Vector computeDiag();
        Vector bmin;
        Vector bmax;
};
#endif