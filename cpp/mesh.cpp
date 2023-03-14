
#include "../h/mesh.h"
#include "meshLoader.cpp"
#include "boundingBox.cpp"
#include "boxNode.cpp"

BoxNode* Mesh::createBoxNode(int startIndex, int endIndex){
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
Mesh::Mesh(string shapeName, Vector albedo, Texture texture) : meshData(shapeName){
    // meshData.readOBJ(file);
    this->root = createBoxNode(0, meshData.indices.size());
    this->albedo = albedo;
    this->texture = texture;
    this->refractionIndex = 1.5;
};
void Mesh::scale(double scale_ratio){
    meshData.scale(scale_ratio);
};
void Mesh::translate(Vector translation){
    meshData.translate(translation);
};

Vector Mesh::relativeTranslate(Vector translation){
    return meshData.relativeTranslate(translation);
};

void Mesh::rotate(Vector center, Vector rotation){
    meshData.rotate(center, rotation);
}

void Mesh::relativeRotate(Vector center, Vector rotation){
    meshData.relativeRotate(center, rotation);
}
void Mesh::updateBoxNode(){
    this->root = createBoxNode(0, meshData.indices.size());
}

Vector Mesh::intersection(Ray ray, float& t, Vector& normal){
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