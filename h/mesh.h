#ifndef MESH_H
#define MESH_H
#include "geometry.h"
#include "meshLoader.h"
#include "boxNode.h"
#include "vector.h"
#include "boundingBox.h"

class Mesh : public Geometry {
    private:
    BoxNode* createBoxNode(int startIndex, int endIndex);
    public:
        Mesh(string shapeName, Vector albedo, Texture texture = Texture::Diffuse);
        void scale(double scale_ratio);
        void translate(Vector translation);
        Vector relativeTranslate(Vector translation);
        void rotate(Vector center, Vector rotation);
        void relativeRotate(Vector center, Vector rotation);
        void updateBoxNode();
        Vector intersection(Ray ray, float& t, Vector& normal);

        TriangleMesh meshData;
        BoxNode* root;
        int textureWidth, textureHeight;
};
#endif