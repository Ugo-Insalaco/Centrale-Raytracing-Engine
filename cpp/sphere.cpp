#include "../h/sphere.h"
Sphere::Sphere(Vector position, Vector albedo, double radius, Texture texture, float refractionIndex){
    this->position = position;
    this->albedo = albedo;
    this->radius = radius;
    this->texture = texture;
    this->refractionIndex = refractionIndex;
}
Vector Sphere::intersection(Ray ray, float& t, Vector& normal){
    double b = 2*dot(ray.direction, ray.origin - this->position);
    double a = (ray.direction).norm2();
    double c = (ray.origin - this->position).norm2() - sqr(radius);
    double delta = sqr(b) - 4*a*c;
    if(delta >=0){
        double sqrtDelta = sqrt(delta);
        double t1 = (-b-sqrtDelta)/(2*a);
        double t2 = (-b+sqrtDelta)/(2*a);
        t = -1;
        if(t1 >=0){
            t = t1;
        }
        if (t2 >=0 && (t2 <= t1 || t1 <= 0)){
            t = t2;
        }

        normal = ray.origin + t * ray.direction - position;
        return this->albedo;
    }
    t = -1;
    return this->albedo;
}