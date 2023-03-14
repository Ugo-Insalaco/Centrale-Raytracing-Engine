
#include "../h/sphere.h"
#include "../h/scene.h"
#include "../h/mesh.h"

Sphere mk_sphere7(Vector(20,25,10), Vector(0.5,0.2,0.5), 5);
Sphere mk_wall1(Vector(0,0,-1000), Vector(0,0,0), 940);
Sphere mk_wall2(Vector(0,1000,0), Vector(0,0,0), 940);
Sphere mk_wall3(Vector(0,0,1000), Vector(0,0,0), 940);
Sphere mk_wall4(Vector(0,-1000,0), Vector(0,0,0), 990);
Sphere mk_wall5(Vector(1000,0,0), Vector(0,0,0), 940);
Sphere mk_wall6(Vector(-1000,0,0), Vector(0,0,0), 940);

Camera mk_camera(Vector(20,11,55), Vector(0.05*PI/8, 1.46*PI/8, 0), PI/48, 1024, 1024);
Scene mk_scene(mk_camera, mk_sphere7);

Mesh mk_rainbow("Rainbow", Vector(0,0,0), Texture::Diffuse);
mk_rainbow.scale(5);
mk_rainbow.translate(Vector(0, 0, 10));
mk_rainbow.updateBoxNode();
mk_scene.addGeometry(&mk_rainbow);

Vector mk_vehicleRotation(0, -0.74*PI/4, 0);

Mesh mk_mantis("Gold Mantis", Vector(0,0,0), Texture::Diffuse);
mk_mantis.scale(0.6);
Vector mk_mantisPosition(-18.2, 9.7, 17);
mk_mantis.translate(mk_mantisPosition);
mk_mantis.rotate(mk_mantisPosition, mk_vehicleRotation);
mk_mantis.updateBoxNode();
mk_scene.addGeometry(&mk_mantis);

Mesh mk_yoshi("Yoshi", Vector(0,0,0), Texture::Diffuse);
mk_yoshi.scale(0.6);
Vector mk_yoshiPosition(-18.2, 10, 17);
mk_yoshi.translate(mk_yoshiPosition);
mk_yoshi.rotate(mk_yoshiPosition, mk_vehicleRotation);
mk_yoshi.updateBoxNode();
mk_scene.addGeometry(&mk_yoshi);

Mesh mk_mushmellow("Mushmellow", Vector(0,0,0), Texture::Diffuse);
mk_mushmellow.scale(0.6);
Vector mk_mushmellowPosition(-16.9, 9.7, 17);
mk_mushmellow.translate(mk_mushmellowPosition);
mk_mushmellow.rotate(mk_mushmellowPosition, mk_vehicleRotation);
mk_mushmellow.updateBoxNode();
mk_scene.addGeometry(&mk_mushmellow);

Mesh mk_waluigi("Waluigi", Vector(0,0,0), Texture::Diffuse);
mk_waluigi.scale(0.6);
Vector mk_waluigiPosition(-16.9, 10, 17);
mk_waluigi.translate(mk_waluigiPosition);
mk_waluigi.rotate(mk_waluigiPosition, mk_vehicleRotation);
mk_waluigi.updateBoxNode();
mk_scene.addGeometry(&mk_waluigi);

Mesh mk_banisher("Banisher", Vector(0,0,0), Texture::Diffuse);
mk_banisher.scale(0.6);
Vector mk_banisherPosition(-17.5, 9.7, 17);
mk_banisher.translate(mk_banisherPosition);
mk_banisher.rotate(mk_banisherPosition, mk_vehicleRotation);
Vector mk_banisherRelTranslate = mk_banisher.relativeTranslate(Vector(0, -1, 13));
mk_banisher.relativeRotate(mk_banisherPosition + mk_banisherRelTranslate, Vector(1.5*PI/4, 0, 0));
mk_banisher.relativeTranslate(Vector(-0.3, 0, 0));
mk_banisher.updateBoxNode();
mk_scene.addGeometry(&mk_banisher);

Mesh mk_daisy("Daisy", Vector(0,0,0), Texture::Diffuse);
mk_daisy.scale(0.6);
Vector mk_daisyPosition(-17.5, 10, 17.5);
mk_daisy.translate(mk_daisyPosition);
mk_daisy.rotate(mk_daisyPosition, mk_vehicleRotation);
Vector mk_daisyRelTranslate = mk_daisy.relativeTranslate(Vector(0, -1, 13));
mk_daisy.relativeRotate(mk_daisyPosition + mk_daisyRelTranslate, Vector(1.5*PI/4, 0, 0));
mk_daisy.relativeTranslate(Vector(0.08, 0, 0));
mk_daisy.updateBoxNode();
mk_scene.addGeometry(&mk_daisy);

mk_scene.addGeometry(&mk_wall1);
mk_scene.addGeometry(&mk_wall2);
mk_scene.addGeometry(&mk_wall3);
mk_scene.addGeometry(&mk_wall4);
mk_scene.addGeometry(&mk_wall5);
mk_scene.addGeometry(&mk_wall6);