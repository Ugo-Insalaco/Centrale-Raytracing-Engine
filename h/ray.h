#ifndef RAY_H
#define RAY_H
#include "vector.h"
class Ray{
    public:
        Ray(Vector& origin, Vector& direction);
        Vector origin;
        Vector direction;
};
#endif