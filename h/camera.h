#ifndef CAMERA_H
#define CAMERA_H
class Camera{
    public: 
        Camera(Vector position, Vector rotation, double fov, int width, int height);
        Vector position;
        Vector rotation;
        double fov;
        int width;
        int height;
        vector<unsigned char> image;
};

#endif