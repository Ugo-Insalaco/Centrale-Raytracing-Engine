#include "../h/camera.h"
Camera::Camera(Vector position, Vector rotation, double fov, int width, int height) : 
    position(position), 
    rotation(rotation),
    fov(fov), 
    width(width), 
    height(height), 
    image(width*height*3, 0){};
