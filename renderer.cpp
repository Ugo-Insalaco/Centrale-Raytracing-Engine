#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

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
bool indirectFlag = false;
bool smoothFlag = false;
int NRays = 1;
int NBounce = 4;
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
        virtual Vector intersection(Ray ray, float& t, Vector& normal) = 0;
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
        Mesh(string shapeName, Vector albedo, Texture texture = Texture::Diffuse) : meshData(shapeName){
            // meshData.readOBJ(file);
            this->root = createBoxNode(0, meshData.indices.size());
            this->albedo = albedo;
            this->texture = texture;
            this->refractionIndex = 1.5;
        };
        void scale(double scale_ratio){
            meshData.scale(scale_ratio);
        };
        void translate(Vector translation){
            meshData.translate(translation);
        };

        Vector relativeTranslate(Vector translation){
            return meshData.relativeTranslate(translation);
        };

        void rotate(Vector center, Vector rotation){
            meshData.rotate(center, rotation);
        }

        void relativeRotate(Vector center, Vector rotation){
            meshData.relativeRotate(center, rotation);
        }
        void updateBoxNode(){
            this->root = createBoxNode(0, meshData.indices.size());
        }

        void inverseNormals(){
            this->meshData.inverseNormals();
        }

        Vector intersection(Ray ray, float& t, Vector& normal){
            t = -1;
            Vector albedo(this->albedo);
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
                                    Material* mat = meshData.materials.at(face.group);
                                    Vector uv = alpha * meshData.uvs[face.uvi] + beta * meshData.uvs[face.uvj] + gamma * meshData.uvs[face.uvk];
                                    uv = Vector(abs(uv[0]), abs(uv[1]), 0);
                                    uv = Vector(uv[0] - (int)uv[0], uv[1] - (int)uv[1], 0);
                                    Vector temp(mat->width - 1, mat->height - 1, 1);
                                    uv = uv * temp;
                                    uv[1] = mat->height - 1 - uv[1];
                                    uv = Vector((int)uv[0], (int)uv[1], 1);
                                    if(mat->width > 0){
                                        albedo = mat->texture[uv[0] + mat->width * uv[1]];
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
            return albedo;
        }

        TriangleMesh meshData;
        BoxNode* root;
        int textureWidth, textureHeight;
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
        Vector intersection(Ray ray, float& t, Vector& normal){
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

        double radius;
};

class Camera{
    public: 
        Camera(Vector position, Vector rotation, double fov, int width, int height) : 
            position(position), 
            rotation(rotation),
            fov(fov), 
            width(width), 
            height(height), 
            image(width*height*3, 0){};

        Vector position;
        Vector rotation;
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
        
        Vector intersect(Ray ray, float& t, int& index, Vector& normal){
            Vector albedo;
            t = -1;
            index = -1;
            Vector normalTemp;
            for (unsigned  sphereIndex = 0; sphereIndex < geometries.size(); sphereIndex++){
                Geometry* sphere = geometries.at(sphereIndex);
                float t2;
                Vector temp_albedo = sphere->intersection(ray, t2, normalTemp);
                if(t2 >= 0 && (t2 <= t || t == -1)){
                    t = t2;
                    index = sphereIndex;
                    normal = normalTemp;
                    albedo = temp_albedo;
                }
            }
            return albedo;
        };

        Vector getColor(const Ray &ray, int numRebond, stack<int>& traversedSpheres){
            Vector noir(0,0,0);
            Vector albedo(0,0,0);
            if(numRebond < 0){
                return noir;
            }
            float t;
            int firstSphere;
            Vector N;
            albedo = intersect(ray, t, firstSphere, N);
            if(t >= 0){
                Geometry* renderedSphere = geometries.at(firstSphere);
                Vector P = ray.origin + t * ray.direction;
                N.normalize();

                switch (renderedSphere->texture)
                {
                case Texture::Diffuse:{
                    Vector directColor = getDirect2(P, N, albedo, lightSphere);
                    Vector indirect(0,0,0);
                    if(indirectFlag){
                        Vector omega_i = random_cos(N);
                        Vector Peps = P + epsilon * N;
                        indirect = getColor(Ray(Peps, omega_i), numRebond -1, traversedSpheres)*albedo;
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
        Vector getDirect2(const Vector &point, const Vector &normal, Vector albedo, Sphere lightSphere){
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
            Vector returnColor = lightI / (4*sqr(PI*lightSphere.radius)) * (shadowed == false) * max(dot(normal, xxprime), 0.) * max(dot(dirLum, xxprime), 0.) / pdf / sqr(d) * albedo;
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

                        rotateVector(ray.direction, Vector(0,0,0), -1*camera.rotation);
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
    Sphere sphere7(Vector(20,25,10), Vector(0.5,0.2,0.5), 5);
    Sphere wall1(Vector(0,0,-1000), Vector(0.6,0.6,0.2), 940);
    Sphere wall2(Vector(0,1000,0), Vector(0.5,0.2,0.1), 940);
    Sphere wall3(Vector(0,0,1000), Vector(0.5,0.3,0.2), 940);
    Sphere wall4(Vector(0,-1000,0), Vector(0.1,0.2,0.5), 990);
    Sphere wall5(Vector(1000,0,0), Vector(0.6,0.2,0.5), 940);
    Sphere wall6(Vector(-1000,0,0), Vector(0.1,0.5,0.5), 940);
    // Sphere wall1(Vector(0,0,-1000), Vector(0,0,0), 940);
    // Sphere wall2(Vector(0,1000,0), Vector(0,0,0), 940);
    // Sphere wall3(Vector(0,0,1000), Vector(0,0,0), 940);
    // Sphere wall4(Vector(0,-1000,0), Vector(0,0,0), 990);
    // Sphere wall5(Vector(1000,0,0), Vector(0,0,0), 940);
    // Sphere wall6(Vector(-1000,0,0), Vector(0,0,0), 940);

    Camera camera(Vector(20,11,55), Vector(0.05*PI/8, 1.46*PI/8, 0), PI/48, 512, 512);
    Scene scene(camera, sphere7);

    Mesh rainbow("Rainbow", Vector(0,0,0), Texture::Diffuse);
    rainbow.scale(5);
    rainbow.translate(Vector(0, 0, 10));
    rainbow.updateBoxNode();
    scene.addGeometry(&rainbow);

    Vector vehicleRotation(0, -0.74*PI/4, 0);

    Mesh mantis("Gold Mantis", Vector(0,0,0), Texture::Diffuse);
    mantis.scale(0.6);
    Vector mantisPosition(-18.2, 9.7, 17);
    mantis.translate(mantisPosition);
    mantis.rotate(mantisPosition, vehicleRotation);
    mantis.updateBoxNode();
    scene.addGeometry(&mantis);

    Mesh yoshi("Yoshi", Vector(0,0,0), Texture::Diffuse);
    yoshi.scale(0.6);
    Vector yoshiPosition(-18.2, 10, 17);
    yoshi.translate(yoshiPosition);
    yoshi.rotate(yoshiPosition, vehicleRotation);
    yoshi.updateBoxNode();
    scene.addGeometry(&yoshi);

    Mesh mushmellow("Mushmellow", Vector(0,0,0), Texture::Diffuse);
    mushmellow.scale(0.6);
    Vector mushmellowPosition(-16.9, 9.7, 17);
    mushmellow.translate(mushmellowPosition);
    mushmellow.rotate(mushmellowPosition, vehicleRotation);
    mushmellow.updateBoxNode();
    scene.addGeometry(&mushmellow);

    Mesh waluigi("Waluigi", Vector(0,0,0), Texture::Diffuse);
    waluigi.scale(0.6);
    Vector waluigiPosition(-16.9, 10, 17);
    waluigi.translate(waluigiPosition);
    waluigi.rotate(waluigiPosition, vehicleRotation);
    waluigi.updateBoxNode();
    scene.addGeometry(&waluigi);

    Mesh banisher("Banisher", Vector(0,0,0), Texture::Diffuse);
    banisher.scale(0.6);
    Vector banisherPosition(-17.5, 9.7, 17);
    banisher.translate(banisherPosition);
    banisher.rotate(banisherPosition, vehicleRotation);
    Vector banisherRelTranslate = banisher.relativeTranslate(Vector(0, -1, 13));
    banisher.relativeRotate(banisherPosition + banisherRelTranslate, Vector(1.5*PI/4, 0, 0));
    banisher.relativeTranslate(Vector(-0.3, 0, 0));
    banisher.updateBoxNode();
    scene.addGeometry(&banisher);

    Mesh daisy("Daisy", Vector(0,0,0), Texture::Diffuse);
    daisy.scale(0.6);
    Vector daisyPosition(-17.5, 10, 17.5);
    daisy.translate(daisyPosition);
    daisy.rotate(daisyPosition, vehicleRotation);
    Vector daisyRelTranslate = daisy.relativeTranslate(Vector(0, -1, 13));
    daisy.relativeRotate(daisyPosition + daisyRelTranslate, Vector(1.5*PI/4, 0, 0));
    daisy.relativeTranslate(Vector(0.08, 0, 0));
    daisy.updateBoxNode();
    scene.addGeometry(&daisy);

    scene.addGeometry(&wall1);
    scene.addGeometry(&wall2);
    scene.addGeometry(&wall3);
    scene.addGeometry(&wall4);
    scene.addGeometry(&wall5);
    scene.addGeometry(&wall6);
    scene.render();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << duration.count() << endl;
    // double temp;
    // cin >> temp;
    return 0;
}