#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <omp.h>
#include <cmath>
#include <iostream>
#include <stack>
#include <random>
using namespace std;

const double M_PI = 3.141592654;

default_random_engine engine[12];
uniform_real_distribution<double> uniform(0.0, 1.0);

static inline double sqr(double x) { return x * x; }

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double& operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector& operator+=(const Vector& v) {
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	double norm2() const {
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}

    double normalize() {
        double n = sqrt(this->norm2());
        coord[0]/= n;
        coord[1]/= n;
        coord[2]/= n;
        return n;
    }

    void print(){
        cout << "x: " << coord[0] << " y: " << coord[1] << " z: " << coord[2] << endl;
    }
	double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, Vector& b){
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(
    Vector& a, double b){
    return Vector(a[0]/b, a[1]/b, a[2]/b);
}

double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}
class Ray{
    public:
        Ray(Vector& origin, Vector& direction): origin(origin), direction(direction){}

        Vector origin;
        Vector direction;
};

enum class Texture {Diffuse, Mirror, Transparent, Light};
class Sphere{
    public:
        Sphere(Vector position, Vector albedo, double radius, Texture texture = Texture::Diffuse, float refractionIndex = 1){
            this->position = position;
            this->albedo = albedo;
            this->radius = radius;
            this->texture = texture;
            this->refractionIndex = refractionIndex;
        }
        float intersection(Ray ray){
            double b = 2*dot(ray.direction, ray.origin - this->position);
            double a = (ray.direction).norm2();
            double c = (ray.origin - this->position).norm2() - sqr(radius);
            double delta = sqr(b) - 4*a*c;
            if(delta >=0){
                double sqrtDelta = sqrt(delta);
                double t1 = (-b-sqrtDelta)/(2*a);
                double t2 = (-b+sqrtDelta)/(2*a);
                double t = -1;
                if(t1 >=0){
                    t = t1;
                }
                if (t2 >=0 && (t2 <= t1 || t1 <= 0)){
                    t = t2;
                }
                return t;
            }
            return -1;
        }

        Vector position;
        Vector albedo;
        double radius;
        Texture texture;
        float refractionIndex;
};

class Camera{
    public: 
        Camera(Vector position, double fov, int width, int height) : 
            position(position), 
            fov(fov), 
            width(width), 
            height(height), 
            image(width*height*3, 0){};

        Vector position;
        double fov;
        int width;
        int height;
        vector<unsigned char> image;
};

class Light{
    public: 
        Light(Vector position, double I){
            this->position = position;
            this->I = I;
        };
        Vector position; 
        double I;
};

class Scene{
    public:
        Scene(Camera camera, Light light): camera(camera), light(light){};

        void addSphere(Sphere sphere){
            if(sphere.texture == Texture::Light) lightSpheres.push_back(spheres.size());
            spheres.push_back(sphere);
        };
        
        void intersect(Ray ray, float& t, int& index){
            t = -1;
            index = -1;
            for (unsigned  sphereIndex = 0; sphereIndex < spheres.size(); sphereIndex++){
                Sphere sphere = spheres.at(sphereIndex);
                float t2 = sphere.intersection(ray);
                if(t2 >= 0 && (t2 <= t || t == -1)){
                    t = t2;
                    index = sphereIndex;
                }
            }
        };

        Vector getColor(const Ray &ray, int numRebond, stack<int>& traversedSpheres){
            Vector noir(0,0,255);
            if(numRebond < 0){
                return noir;
            }
            float t;
            int firstSphere;
            intersect(ray, t, firstSphere);
            if(t >= 0){
                Sphere renderedSphere = spheres.at(firstSphere);
                Vector P = ray.origin + t * ray.direction;
                Vector N = (P - renderedSphere.position);
                N.normalize();

                switch (renderedSphere.texture)
                {
                case Texture::Light:{
                    return noir;
                }
                case Texture::Diffuse:{
                    // Vector directColor = getDirect(P, N, renderedSphere);
                    Vector directColor = getDirect2(P, N, renderedSphere, spheres[lightSpheres[0]]);
                    Vector omega_i = random_cos(N);
                    Vector Peps = P + epsilon * N;
                    Vector indirect = getColor(Ray(Peps, omega_i), numRebond -1, traversedSpheres)*renderedSphere.albedo;
                    return directColor + indirect;
                    break;
                }
                case Texture::Mirror:{
                    return getBounce(P, N, ray, numRebond, traversedSpheres);
                    break;
                }
                case Texture::Transparent: {
                    float n1, n2;
                    bool in;
                    if(traversedSpheres.empty()){ // On rentre dans la première sphère
                        n1 = ambiantIndex;
                        n2 = renderedSphere.refractionIndex;
                        in = true;
                    }
                    else{
                        n1 = spheres[traversedSpheres.top()].refractionIndex;
                        if(traversedSpheres.top() == firstSphere){ // Si on sort de la sphère
                            in = false;
                            N = (-1)*N;
                            traversedSpheres.pop();
                            if(traversedSpheres.empty()){ // Si on sort de la dernière sphère de la pile
                                n2 = ambiantIndex;
                            }
                            else{ // S'il reste des sphères dans la pile
                                n2 = spheres[traversedSpheres.top()].refractionIndex;
                            }
                        }
                        else{ // Si on rentre dans la sphère
                            in = true;
                            n2 = renderedSphere.refractionIndex;
                        }
                    }

                    float dotN = dot(N, ray.direction);
                    float k0 = pow(n1-n2, 2)/pow(n1+n2, 2);
                    float R = k0 + (1 - k0)*pow((1-abs(dotN)), 5);
                    int thread_id = omp_get_thread_num();
                    double r = uniform(engine[thread_id]);
                    float decision = 1-sqr(n1/n2)*(1-sqr(dotN));

                    if(decision > 0 && r > R){
                        if(in){
                            traversedSpheres.push(firstSphere);
                        }
                        Vector wtT = (n1/n2)*(ray.direction+(-dotN)*N);
                        Vector wtN = -sqrt(decision)*N;
                        Vector wt = wtT + wtN;
                        Vector Peps = P + (-1)*epsilon * N;
                        Ray newRay(Peps, wt);
                        return getColor(newRay, numRebond - 1, traversedSpheres);
                    }
                    return getBounce(P, N, ray, numRebond, traversedSpheres);
                    break;
                }
                default: 
                    return noir;
                    break;
                }
            };
            return noir;
        }
        Vector random_cos(const Vector& N){
            int thread_id = omp_get_thread_num();
            double r1 = uniform(engine[thread_id]);
            double r2 = uniform(engine[thread_id]);
            double r = sqrt(1 - r2);
            double x = r*cos(2. * M_PI * r1);
            double y = r*sin(2. * M_PI * r1);
            double z = sqrt(r2);

            Vector T1;
            if((abs(N[0]) <= abs(N[1])) && (abs(N[0]) <= abs(N[2]))){
                T1 = Vector(0, -N[2], N[1]);
            }
            if((abs(N[1]) <= abs(N[2])) && (abs(N[1]) <= abs(N[0]))){
                T1 = Vector(-N[2], 0, N[0]);
            }
            else{
                T1 = Vector(-N[1], N[0], 0);
            }
            T1.normalize();
            Vector T2 = cross(N, T1);
            T2.normalize();
            return z*N+x*T1+y*T2;
        }
        Vector getBounce(const Vector &point, const Vector &normal, const Ray &ray, int numRebond, stack <int> &traversedSpheres){
            Vector newDirection = ray.direction - 2*dot(ray.direction, normal)*normal;
            Vector Peps = point + epsilon * normal;
            Ray newRay(Peps, newDirection);
            return getColor(newRay, numRebond - 1, traversedSpheres);
        };
        Vector getDirect2(const Vector &point, const Vector &normal, Sphere &renderedSphere, Sphere &lightSphere){
            Vector dirLum = lightSphere.position - point;
            dirLum.normalize();
            Vector dir_xprime = random_cos(-1 * dirLum);
            Vector xprime = dir_xprime*lightSphere.radius + lightSphere.position;
            Vector xxprime = xprime - point;
            double d = xxprime.norm2();
            xxprime.normalize();

            Vector Peps = point + epsilon * normal;
            Ray lightRay(Peps, xxprime);
            bool shadowed = false;

            float otherSphereIntersection;
            int i;
            intersect(lightRay, otherSphereIntersection, i);
            if(otherSphereIntersection > 0 && sqr(otherSphereIntersection)*lightRay.direction.norm2() < d && i != lightSpheres[0]){
                shadowed = true;
            }
            double pdf = -dot(dirLum, dir_xprime)/(M_PI * sqr(lightSphere.radius));
            // cout << (shadowed == false) << endl;
            Vector returnColor = light.I / (4*sqr(M_PI*lightSphere.radius)) * (shadowed == false) * max(dot(normal, xxprime), 0.) * max(dot(dirLum, xxprime), 0.) / pdf / sqr(d) * renderedSphere.albedo;
            return returnColor;
        }
        Vector getDirect(const Vector &point, const Vector &normal, Sphere &renderedSphere){
            Vector wi = (light.position - point);
            double d = wi.normalize();

            Vector Peps = point + epsilon * normal;
            Ray lightRay(Peps, wi);
            bool shadowed = false;

            float otherSphereIntersection;
            int i;
            intersect(lightRay, otherSphereIntersection, i);
            if(otherSphereIntersection > 0 && sqr(otherSphereIntersection)*lightRay.direction.norm2() < sqr(d)){
                shadowed = true;
            }

            double const c = light.I/(4*sqr(M_PI)*sqr(d))*dot(normal, wi) * (shadowed == 0) * (dot(normal, lightRay.direction) >= 0);
            Vector colors = c/M_PI * renderedSphere.albedo;
            return colors;
        }

        void render(){
            int pixel_done = 0;
            #pragma omp parallel for
            for (int i = 0; i < camera.height; i++) {
                cout << (float)(pixel_done)/(camera.height*camera.width) * 100<<"%" << endl;
                for (int j = 0; j < camera.width; j++) {
                    pixel_done++;
                    int NRays = 500;
                    int NBounce = 5;
                    Vector finalColor;
                    double dfocus = 50;
                    for (int rayId = 0; rayId < NRays; rayId++){   
                        int thread_id = omp_get_thread_num();
                        double r1 = uniform(engine[thread_id]);
                        double r2 = uniform(engine[thread_id]);
                        double r = sqrt(-2*log(r1));
                        double gx = r*cos(2*M_PI*r2)*0.7;
                        double gy = r*sin(2*M_PI*r2)*0.7;

                        r1 = uniform(engine[thread_id]);
                        r2 = uniform(engine[thread_id]);
                        r = sqrt(-2*log(r1));
                        double cx = r*cos(2*M_PI*r2)*0.7;
                        double cy = r*sin(2*M_PI*r2)*0.7;

                        Vector pixel(j-camera.width/2 + 0.5 + gx, -i +camera.height/2 - 0.5 + gy, -camera.width/(2*tan(camera.fov/2)));
                        pixel.normalize();
                        Vector Pp = camera.position + dfocus * pixel;
                        Vector Op = camera.position + Vector(cx, cy, 0);
                        Vector newDir = (Pp + (-1) * Op);
                        newDir.normalize();
                        Ray ray(Op,newDir);
                        stack<int> traversed;
                        Vector color = getColor(ray,NBounce , traversed);
                        finalColor += color;
                    }
                    finalColor = finalColor / NRays;
                    camera.image[(i*camera.width + j) * 3 + 0] = min(int(pow(finalColor[0], 1/2.2)), 255);   // RED
                    camera.image[(i*camera.width + j) * 3 + 1] = min(int(pow(finalColor[1], 1/2.2)), 255);  // GREEN
                    camera.image[(i*camera.width + j) * 3 + 2] = min(int(pow(finalColor[2], 1/2.2)), 255);  // BLUE
                    
                }
            }
            stbi_write_png("image.png", camera.width, camera.height, 3, &(camera.image[0]), 0);
        };
    private:
        vector<Sphere> spheres;
        vector<int> lightSpheres;
        Camera camera;
        Light light;
        const float epsilon = 0.01;
        const float ambiantIndex = 1;
};

int main() {
    Sphere sphere(Vector(0,6,20), Vector(0.5,0.2,0.5), 6);
    Sphere sphere2(Vector(-15,12,10), Vector(0.5,0.2,0.5), 12, Texture::Mirror);
    Sphere sphere3(Vector(0,20,0), Vector(0.5,0.2,0.5), 5, Texture::Transparent, 1.3);
    // Sphere sphere6(Vector(0,18,0), Vector(0.5,0.2,0.5), 6);
    Sphere sphere7(Vector(0,24,50), Vector(0.5,0.2,0.5), 2, Texture::Light);
    // Sphere sphere4(Vector(-6,-3,0), Vector(0.5,0.2,0.5), 5);
    // Sphere sphere5(Vector(8,-3,0), Vector(0.5,0.2,0.5), 6);
    Sphere wall1(Vector(0,0,-1000), Vector(0.6,0.6,0.2), 940);
    Sphere wall2(Vector(0,1000,0), Vector(0.5,0.2,0.1), 940);
    Sphere wall3(Vector(0,0,1000), Vector(0.5,0.3,0.2), 940);
    Sphere wall4(Vector(0,-1000,0), Vector(0.1,0.2,0.5), 990);
    Sphere wall5(Vector(1000,0,0), Vector(0.6,0.2,0.5), 940);
    Sphere wall6(Vector(-1000,0,0), Vector(0.1,0.5,0.5), 940);
    Camera camera(Vector(0,10,55), M_PI/3, 1024, 1024);
    Light light(Vector(-10, 20, 40), 2E12);

    Scene scene(camera, light);
    scene.addSphere(sphere);
    scene.addSphere(sphere2);
    scene.addSphere(sphere3);
    // scene.addSphere(sphere4);
    // scene.addSphere(sphere5);
    // scene.addSphere(sphere6);
    scene.addSphere(sphere7);
    scene.addSphere(wall1);
    scene.addSphere(wall2);
    scene.addSphere(wall3);
    scene.addSphere(wall4);
    scene.addSphere(wall5);
    scene.addSphere(wall6);
    scene.render();

    return 0;
}