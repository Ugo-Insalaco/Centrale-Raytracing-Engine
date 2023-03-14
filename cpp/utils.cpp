#include "../h/boundingBox.h"
#include "../h/meshLoader.h"
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