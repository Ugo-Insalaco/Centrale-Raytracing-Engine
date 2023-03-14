#include "../h/boundingBox.h"
BoundingBox::BoundingBox(){};
BoundingBox::BoundingBox(Vector& bmin, Vector& bmax) : bmin(bmin), bmax(bmax) {};
bool BoundingBox::intersection(Ray ray){
    Vector N, P;
    double dotu;
    N = Vector(1, 0, 0);
    dotu = 1/dot(N, ray.direction);
    P = Vector(bmin[0], 0, 0);
    double tminX = dot(P - ray.origin, N)*dotu;
    P = Vector(bmax[0], 0, 0);
    double tmaxX = dot(P - ray.origin, N)*dotu;
    if(tmaxX < tminX) {
        swap(tmaxX, tminX);
    }

    N = Vector(0, 1, 0);
    dotu = 1/dot(N, ray.direction);
    P = Vector(0, bmin[1], 0);
    double tminY = dot(P - ray.origin, N)*dotu;
    P = Vector(0, bmax[1], 0);
    double tmaxY = dot(P - ray.origin, N)*dotu;
    if(tmaxY < tminY) {
        swap(tmaxY, tminY);
    }

    N = Vector(0, 0, 1);
    dotu = 1/dot(N, ray.direction);
    P = Vector(0, 0, bmin[2]);
    double tminZ = dot(P - ray.origin, N)*dotu;
    P = Vector(0, 0, bmax[2]);
    double tmaxZ = dot(P - ray.origin, N)*dotu;
    if(tmaxZ < tminZ) {
        swap(tmaxZ, tminZ);
    }
    double t0 = max(tminX, max(tminY, tminZ)), t1 = min(tmaxX, min(tmaxY, tmaxZ));
    return t0 < t1;
};
Vector BoundingBox::computeDiag(){
    return bmax - bmin;
};