#ifndef SPHERE_H
#define SPHERE_H
class Sphere : public Geometry {
    public:
        Sphere(Vector position, Vector albedo, double radius, Texture texture = Texture::Diffuse, float refractionIndex = 1);
        Vector intersection(Ray ray, float& t, Vector& normal);
        double radius;
};
#endif