#include "../h/sphere.h"
#include "../h/scene.h"
#include "../h/mesh.h"

Sphere c_sphere7(Vector(20,25,10), Vector(0.5,0.2,0.5), 5);
Sphere c_wall1(Vector(0,0,-1000), Vector(0.7,0.4,0.4), 940);
Sphere c_wall2(Vector(0,1000,0), Vector(0,0.1,0), 940);
Sphere c_wall3(Vector(0,0,1000), Vector(0,0,0.2), 940);
Sphere c_wall4(Vector(0,-1000,0), Vector(0,0.1,0.6), 990);
Sphere c_wall5(Vector(1000,0,0), Vector(0.3,0,0.4), 940);
Sphere c_wall6(Vector(-1000,0,0), Vector(0.6,0.8, 0.2), 940);

Camera c_camera(Vector(0,0,55), Vector(0,0, 0), PI/3, 256, 256);
Scene c_scene(c_camera, c_sphere7);

Mesh c_daisy("daisy", Vector(0,0,0), Texture::Diffuse);
c_daisy.scale(10);
c_daisy.translate(Vector(0, -10, 0));
c_daisy.updateBoxNode();

Mesh c_toad("yoshi", Vector(0,0,0), Texture::Diffuse);
c_toad.scale(0.1);
c_toad.translate(Vector(15, -10, 0));
c_toad.updateBoxNode();


c_scene.addGeometry(&c_wall1);
c_scene.addGeometry(&c_wall2);
c_scene.addGeometry(&c_wall3);
c_scene.addGeometry(&c_wall4);
c_scene.addGeometry(&c_wall5);
c_scene.addGeometry(&c_wall6);
c_scene.addGeometry(&c_daisy);
c_scene.addGeometry(&c_toad);