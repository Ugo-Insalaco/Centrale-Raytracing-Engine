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
#include <list>
#include <chrono>
#include "meshLoader.cpp"
using namespace std;

default_random_engine engine[4];
uniform_real_distribution<double> uniform(0.0, 1.0);

float PI = 3.141592654;
bool dephtFlag = false;
bool aliasFlag = true;
bool fresnelFlag = true;
bool indirectFlag = true;
bool smoothFlag = false;
int NRays = 10;
int NBounce = 6;
double dfocus = 5;
double lightI = 1E13;

class Ray{
    public:
        Ray(Vector& origin, Vector& direction): origin(origin), direction(direction){}

        Vector origin;
        Vector direction;
};

enum class Texture {Diffuse, Mirror, Transparent};

class BoundingBox{
    public :
        BoundingBox(){};
        BoundingBox(Vector& bmin, Vector& bmax) : bmin(bmin), bmax(bmax) {};
        bool intersection(Ray ray){
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
        Vector computeDiag(){
            return bmax - bmin;
        }
        Vector bmin;
        Vector bmax;
};
class Geometry{
    public:
        virtual void intersection(Ray ray, float& t, Vector& normal) = 0;
        Vector position;
        Vector albedo;
        Texture texture;
        float refractionIndex;
};

BoundingBox computeBoundingBox(TriangleMesh mesh, int startIndex, int endIndex){
    Vector bmin = mesh.vertices[mesh.indices[0].vtxi];
    Vector bmax = mesh.vertices[mesh.indices[0].vtxi];
    for(int faceIndex = startIndex; faceIndex < endIndex; faceIndex++){
        TriangleIndices face = mesh.indices[faceIndex];
        Vector A = mesh.vertices[face.vtxi];
        Vector B = mesh.vertices[face.vtxj];
        Vector C = mesh.vertices[face.vtxk];

        double minX = min(A[0], min(B[0], C[0]));
        double minY = min(A[1], min(B[1], C[1]));
        double minZ = min(A[2], min(B[2], C[2]));
        bmin[0] = min(bmin[0], minX);
        bmin[1] = min(bmin[1], minY);
        bmin[2] = min(bmin[2], minZ);

        double maxX = max(A[0], max(B[0], C[0]));
        double maxY = max(A[1], max(B[1], C[1]));
        double maxZ = max(A[2], max(B[2], C[2]));
        bmax[0] = max(bmax[0], maxX);
        bmax[1] = max(bmax[1], maxY);
        bmax[2] = max(bmax[2], maxZ);
    }
    return BoundingBox(bmin, bmax);
};

class BoxNode {
    public:
        BoxNode(BoundingBox boundingBox, int startIndex, int endIndex) : boundingBox(boundingBox), startIndex(startIndex), endIndex(endIndex){};
        void print(){
            cout << "BoxNode: startindex = " << startIndex << " endIndex = " << endIndex << endl;
            if(leftChild && rightChild){
                leftChild->print();
                rightChild->print();
            }
        };
        void printBoundaries(){
            cout << "BoxNode boundaries: startindex = " << startIndex << " endIndex = " << endIndex << endl;
        }
        BoundingBox boundingBox;
        BoxNode* leftChild {};
        BoxNode* rightChild {};
        int startIndex, endIndex;
};

class Mesh : public Geometry {
    private:
    BoxNode* createBoxNode(int startIndex, int endIndex){
        cout << "startIndex: " << startIndex << " endIndex: " << endIndex << endl;
        BoundingBox boundingBox = computeBoundingBox(meshData, startIndex, endIndex);
        BoxNode* node = new BoxNode(boundingBox, startIndex, endIndex);
        Vector diag = boundingBox.computeDiag();
        Vector middleDiag = boundingBox.bmin + 0.5 * diag;
        int longestAxis = diag.getLongestAxis();
        int pivotIndex = startIndex;
        for(int faceIndex = startIndex; faceIndex < endIndex; faceIndex++){
            TriangleIndices face = meshData.indices[faceIndex];
            Vector A = meshData.vertices[face.vtxi];
            Vector B = meshData.vertices[face.vtxj];
            Vector C = meshData.vertices[face.vtxk];
            Vector faceBarycenter = computeBarycenter(A, B, C);
            if(faceBarycenter[longestAxis] < middleDiag[longestAxis]){
                swap(meshData.indices[faceIndex], meshData.indices[pivotIndex]);
                pivotIndex++;
            }
        }
        if( pivotIndex<=startIndex || pivotIndex>=endIndex-1 || endIndex-startIndex<5 ){
            return node;
        }
        node->leftChild = createBoxNode(startIndex, pivotIndex);
        node->rightChild = createBoxNode(pivotIndex, endIndex);
        return node;
    };
    public:
        Mesh(const char* file, Vector albedo, Texture texture = Texture::Diffuse){
            meshData.readOBJ(file);
            this->root = createBoxNode(0, meshData.indices.size());
            this->albedo = albedo;
            this->texture = texture;
            this->refractionIndex = 1.5;
            root->print();
        };
        void scale(double scale_ratio){
            meshData.scale(scale_ratio);
        };
        void translate(Vector translation){
            meshData.translate(translation);
        };

        void rotate(Vector center, Vector rotation){
            meshData.rotate(center, rotation);
        }
        void updateBoxNode(){
            this->root = createBoxNode(0, meshData.indices.size());
        }
        void intersection(Ray ray, float& t, Vector& normal){
            t = -1;
            if(root->boundingBox.intersection(ray)){
                list<BoxNode*> nodesToVisit;
                nodesToVisit.push_front(root);
                double bestt = numeric_limits<double>::max();
                while(!nodesToVisit.empty()){
                    BoxNode* currentNode = nodesToVisit.back();
                    nodesToVisit.pop_back();
                    if(currentNode->leftChild){
                        if(currentNode->leftChild->boundingBox.intersection(ray)){
                            nodesToVisit.push_back(currentNode->leftChild);
                        }
                        if(currentNode->rightChild->boundingBox.intersection(ray)){
                            nodesToVisit.push_back(currentNode->rightChild);
                        }
                    }
                    else{
                        for(int faceIndex = currentNode->startIndex; faceIndex < currentNode->endIndex; faceIndex++){
                            TriangleIndices face = meshData.indices[faceIndex];
                            Vector A = meshData.vertices[face.vtxi];
                            Vector B = meshData.vertices[face.vtxj];
                            Vector C = meshData.vertices[face.vtxk];
                            Vector e1 = B - A;
                            Vector e2 = C - A;
                            Vector N = cross(e1, e2);
                            double den = dot(ray.direction, N);
                            double beta = dot(e2, cross(A - ray.origin, ray.direction))/den;
                            double gamma = -dot(e1, cross(A - ray.origin, ray.direction))/den;
                            double alpha = 1 - beta - gamma;
                            double ttemp = dot(A - ray.origin, N)/den;
                            if(alpha > 0 && alpha < 1 && beta > 0 && beta < 1 && gamma > 0 && gamma < 1 && ttemp > 0){
                                if(ttemp <= bestt){
                                    if(smoothFlag){
                                         Vector NA = meshData.normals[face.vtxi];
                                        Vector NB = meshData.normals[face.vtxj];
                                        Vector NC = meshData.normals[face.vtxk];
                                        N = alpha*NA+beta*NB+gamma*NC;
                                        N.normalize();
                                    }
                                    bestt = ttemp;
                                    normal = N;
                                }
                            }
                        }
                    }
                }
                t = bestt;
            }
        }

        TriangleMesh meshData;
        BoxNode* root;
};

class Sphere : public Geometry {
    public:
        Sphere(Vector position, Vector albedo, double radius, Texture texture = Texture::Diffuse, float refractionIndex = 1){
            this->position = position;
            this->albedo = albedo;
            this->radius = radius;
            this->texture = texture;
            this->refractionIndex = refractionIndex;
        }
        void intersection(Ray ray, float& t, Vector& normal){
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
                return;
            }
            t = -1;
        }

        double radius;
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


class Scene{
    public:
        Scene(Camera camera, Sphere lightSphere): lightSphere(lightSphere), camera(camera){};

        void addGeometry(Geometry* geometry){
            geometries.push_back(geometry);
        };
        
        void intersect(Ray ray, float& t, int& index, Vector& normal){
            t = -1;
            index = -1;
            Vector normalTemp;
            for (unsigned  sphereIndex = 0; sphereIndex < geometries.size(); sphereIndex++){
                Geometry* sphere = geometries.at(sphereIndex);
                float t2;
                sphere->intersection(ray, t2, normalTemp);
                if(t2 >= 0 && (t2 <= t || t == -1)){
                    t = t2;
                    index = sphereIndex;
                    normal = normalTemp;
                }
            }
        };

        Vector getColor(const Ray &ray, int numRebond, stack<int>& traversedSpheres){
            Vector noir(0,0,0);
            if(numRebond < 0){
                return noir;
            }
            float t;
            int firstSphere;
            Vector N;
            intersect(ray, t, firstSphere, N);
            if(t >= 0){
                Geometry* renderedSphere = geometries.at(firstSphere);
                Vector P = ray.origin + t * ray.direction;
                N.normalize();

                switch (renderedSphere->texture)
                {
                case Texture::Diffuse:{
                    Vector directColor = getDirect2(P, N, renderedSphere, lightSphere);
                    Vector indirect(0,0,0);
                    if(indirectFlag){
                        Vector omega_i = random_cos(N);
                        Vector Peps = P + epsilon * N;
                        indirect = getColor(Ray(Peps, omega_i), numRebond -1, traversedSpheres)*renderedSphere->albedo;
                    }
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
                        n2 = renderedSphere->refractionIndex;
                        in = true;
                    }
                    else{
                        n1 = geometries[traversedSpheres.top()]->refractionIndex;
                        if(traversedSpheres.top() == firstSphere){ // Si on sort de la sphère
                            in = false;
                            N = (-1)*N;
                            traversedSpheres.pop();
                            if(traversedSpheres.empty()){ // Si on sort de la dernière sphère de la pile
                                n2 = ambiantIndex;
                            }
                            else{ // S'il reste des sphères dans la pile
                                n2 = geometries[traversedSpheres.top()]->refractionIndex;
                            }
                        }
                        else{ // Si on rentre dans la sphère
                            in = true;
                            n2 = renderedSphere->refractionIndex;
                        }
                    }

                    float dotN = dot(N, ray.direction);
                    float R;
                    if(fresnelFlag){
                        float k0 = pow(n1-n2, 2)/pow(n1+n2, 2);
                        R = k0 + (1 - k0)*pow((1-abs(dotN)), 5);
                    }
                    else{
                        R = 0;
                    }

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
            double x = r*cos(2. * PI * r1);
            double y = r*sin(2. * PI * r1);
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
        Vector getDirect2(const Vector &point, const Vector &normal, Geometry* renderedGeometry, Sphere lightSphere){
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
            Vector N;
            intersect(lightRay, otherSphereIntersection, i, N);
            if(otherSphereIntersection > 0 && sqr(otherSphereIntersection)*lightRay.direction.norm2() < d){
                shadowed = true;
            }
            double pdf = -dot(dirLum, dir_xprime)/(PI * sqr(lightSphere.radius));
            Vector returnColor = lightI / (4*sqr(PI*lightSphere.radius)) * (shadowed == false) * max(dot(normal, xxprime), 0.) * max(dot(dirLum, xxprime), 0.) / pdf / sqr(d) * renderedGeometry->albedo;
            return returnColor;
        }

        void render(){
            int pixel_done = 0;
            #pragma omp parallel for
            for (int i = 0; i < camera.height; i++) {
                cout << (float)(pixel_done)/(camera.height*camera.width) * 100<<"%" << endl;
                for (int j = 0; j < camera.width; j++) {
                    pixel_done++;
                    Vector finalColor;
                    for (int rayId = 0; rayId < NRays; rayId++){   
                        int thread_id = omp_get_thread_num();
                        Vector pixel;
                        if(aliasFlag){
                            double r1 = uniform(engine[thread_id]);
                            double r2 = uniform(engine[thread_id]);
                            double r = sqrt(-2*log(r1));
                            double gx = r*cos(2*PI*r2)*0.7;
                            double gy = r*sin(2*PI*r2)*0.7;
                            pixel = Vector(j-camera.width/2 + 0.5 + gx, -i +camera.height/2 - 0.5 + gy, -camera.width/(2*tan(camera.fov/2)));
                        }
                        else{
                            pixel = Vector(j-camera.width/2 + 0.5, -i +camera.height/2 - 0.5, -camera.width/(2*tan(camera.fov/2)));
                        }

                        pixel.normalize();
                        
                        Ray ray(camera.position, pixel);
                        if(dephtFlag){
                            double r1 = uniform(engine[thread_id]);
                            double r2 = uniform(engine[thread_id]);
                            double r = sqrt(-2*log(r1));
                            double cx = r*cos(2*PI*r2)*0.7;
                            double cy = r*sin(2*PI*r2)*0.7;

                            Vector Pp = camera.position + dfocus * pixel;
                            Vector Op = camera.position + Vector(cx, cy, 0);
                            Vector newDir = (Pp + (-1) * Op);
                            newDir.normalize();
                            ray = Ray(Op,newDir);
                        }

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
        vector<Geometry*> geometries;
        Sphere lightSphere;
        Camera camera;
        // Light light;
        const float epsilon = 0.01;
        const float ambiantIndex = 1;
};

int main() {
    auto start = chrono::high_resolution_clock::now();
    Sphere sphere(Vector(0,6,20), Vector(0.5,0.2,0.5), 6);
    Sphere sphere2(Vector(-15,12,10), Vector(0.5,0.2,0.5), 12, Texture::Mirror);
    Sphere sphere3(Vector(0,20,0), Vector(0.5,0.2,0.5), 5, Texture::Transparent, 1.3);
    Sphere sphere7(Vector(0,30,-5), Vector(0.5,0.2,0.5), 5);
    Sphere wall1(Vector(0,0,-1000), Vector(0.6,0.6,0.2), 940);
    Sphere wall2(Vector(0,1000,0), Vector(0.5,0.2,0.1), 940);
    Sphere wall3(Vector(0,0,1000), Vector(0.5,0.3,0.2), 940);
    Sphere wall4(Vector(0,-1000,0), Vector(0.1,0.2,0.5), 990);
    Sphere wall5(Vector(1000,0,0), Vector(0.6,0.2,0.5), 940);
    Sphere wall6(Vector(-1000,0,0), Vector(0.1,0.5,0.5), 940);
    Camera camera(Vector(0,10,55), PI/3, 512, 512);
    Mesh mesh("./cat/cat.obj", Vector(0,0.5,0.3), Texture::Mirror);
    mesh.scale(0.5);
    mesh.rotate(Vector(0,0,0), Vector(0, PI/2, 0));
    mesh.translate(Vector(0, 0, 40));
    mesh.updateBoxNode();
    Scene scene(camera, sphere7);
    scene.addGeometry(&wall1);
    scene.addGeometry(&wall2);
    scene.addGeometry(&wall3);
    scene.addGeometry(&wall4);
    scene.addGeometry(&wall5);
    scene.addGeometry(&wall6);
    scene.addGeometry(&mesh);
    scene.render();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << duration.count() << endl;
    double temp;
    cin >> temp;
    return 0;
}