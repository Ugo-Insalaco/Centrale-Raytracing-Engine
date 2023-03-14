#include "../h/scene.h"
default_random_engine engine[4];
uniform_real_distribution<double> uniform(0.0, 1.0);

Scene::Scene(): lightSphere(Sphere(Vector(0,0,0), Vector(0,0,0), 1.)), camera(Camera(Vector(0,0,0), Vector(0,0,0), PI/2, 512, 512)){};
Scene::Scene(Camera camera, Sphere lightSphere): lightSphere(lightSphere), camera(camera){};

void Scene::addGeometry(Geometry* geometry){
    geometries.push_back(geometry);
};

Vector Scene::intersect(Ray ray, float& t, int& index, Vector& normal){
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

Vector Scene::getColor(const Ray &ray, int numRebond, stack<int>& traversedSpheres){
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
Vector Scene::random_cos(const Vector& N){
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
Vector Scene::getBounce(const Vector &point, const Vector &normal, const Ray &ray, int numRebond, stack <int> &traversedSpheres){
    Vector newDirection = ray.direction - 2*dot(ray.direction, normal)*normal;
    Vector Peps = point + epsilon * normal;
    Ray newRay(Peps, newDirection);
    return getColor(newRay, numRebond - 1, traversedSpheres);
};
Vector Scene::getDirect2(const Vector &point, const Vector &normal, Vector albedo, Sphere lightSphere){
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

void Scene::render(){
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