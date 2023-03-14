#ifndef SCENE_H
#define SCENE_H
class Scene{
    public:
        Scene();
        Scene(Camera camera, Sphere lightSphere);
        void addGeometry(Geometry* geometry);
        Vector intersect(Ray ray, float& t, int& index, Vector& normal);
        Vector getColor(const Ray &ray, int numRebond, stack<int>& traversedSpheres);
        Vector random_cos(const Vector& N);
        Vector getBounce(const Vector &point, const Vector &normal, const Ray &ray, int numRebond, stack <int> &traversedSpheres);
        Vector getDirect2(const Vector &point, const Vector &normal, Vector albedo, Sphere lightSphere);
        void render();
        vector<Geometry*> geometries;
        Sphere lightSphere;
        Camera camera;
        const float epsilon = 0.01;
        const float ambiantIndex = 1;
};
#endif