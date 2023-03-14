#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "../h/stb_image.h"
#endif

#include <string>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>
#include <cmath>
#include "../h/meshLoader.h"
#include "triangleIndices.cpp"
#include <cstring>
using namespace std;
 
TriangleMesh::TriangleMesh(string shapeName) {
    this->rotation = Vector(0,0,0);
    this->path = "./Models/"+shapeName;
    this->fileName = shapeName;
    string mtlFile = this->path+"/"+this->fileName+".mtl";
    string objFile = this->path+"/"+this->fileName+".obj";
    readMTL(mtlFile);
    readOBJ(objFile);
    loadTextures(path);
};
void TriangleMesh::readMTL(string mtl) {
    char grp[255];
    char mapFile[255];
    FILE* f;
    f = fopen(mtl.c_str(), "r");
    string curGroup = "";
    Material* mat = new Material();

    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'n' && line[1] == 'e' && line[2] == 'w') {
            if(curGroup != ""){
                materials.insert(pair<string, Material*>(curGroup, mat));
            }
            mat = new Material();
            sscanf(line, "newmtl %[^\n]\n", grp);
            curGroup = grp;
        }

        if (line[0] == 'K' && line[1] == 'd') {
            Vector vec;
            sscanf(line, "Kd %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            mat->Kd = vec;
        }

        if (line[0] == 'm' && line[1] == 'a' && line[2] == 'p') {
            sscanf(line, "map_Kd %[^\n]\n", mapFile);
            mat->KdMap = mapFile;
        }

    }
    materials.insert(pair<string, Material*>(curGroup, mat));
    fclose(f);
}
void TriangleMesh::readOBJ(string obj) {

    // char matfile[255];
    char grp[255];

    FILE* f;
    f = fopen(obj.c_str(), "r");
    string curGroup = "";
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup = grp;
        }

        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                col[0] = min(1., max(0., col[0]));
                col[1] = min(1., max(0., col[1]));
                col[2] = min(1., max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);

            } else {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char* consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                indices.push_back(t);
            } else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true) {
                if (consumedline[0] == '\n') break;
                if (consumedline[0] == '\0') break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                    if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                    if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                    if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                    if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                    if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                    if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                } else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    } else {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            } else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }

        }

    }
    fclose(f);
}

void TriangleMesh::scale(double scale_ratio){
    for(unsigned i=0; i<vertices.size(); i++){
        vertices[i] = scale_ratio*vertices[i];
    }
}

void TriangleMesh::translate(Vector translation){
    for(unsigned i=0; i<vertices.size(); i++){
        vertices[i] = vertices[i] + translation;
    }
}

Vector TriangleMesh::relativeTranslate(Vector translation){
    Vector newTranslation = translation;
    rotateVector(newTranslation, Vector(0,0,0), rotation);
    translate(newTranslation);
    return newTranslation;
}

void TriangleMesh::rotate(Vector center, Vector rotation){
    float alpha = rotation[0];
    float beta = rotation[1];
    float gamma = rotation[2];
    this->rotation += rotation;
    translate(-1*center);
    for(unsigned i=0; i<vertices.size(); i++){
        Vector v = vertices[i];
        vertices[i] = Vector(v[0], cos(alpha)*v[1]-sin(alpha)*v[2], sin(alpha)*v[1]+cos(alpha)*v[2]); 
    }
    for(unsigned i=0; i<vertices.size(); i++){
        Vector v = vertices[i];
        vertices[i] = Vector(cos(beta)*v[0]-sin(beta)*v[2], v[1], sin(beta)*v[0]+cos(beta)*v[2]); 
    }
    for(unsigned i=0; i<vertices.size(); i++){
        Vector v = vertices[i];
        vertices[i] = Vector(cos(gamma)*v[0]-sin(gamma)*v[1], sin(gamma)*v[0]+cos(gamma)*v[1], v[2]); 
    }
    translate(center);
}

void TriangleMesh::relativeRotate(Vector center, Vector rotation){
    Vector rot = this->rotation;
    rotate(center, -1*this->rotation);
    rotate(center, rotation);
    rotate(center, rot);
}

void TriangleMesh::loadTextures(string path){
    map<string, Material*>::iterator it;
    for (it = materials.begin(); it != materials.end(); it++)
    {
        if(it->second->KdMap != ""){
            int width, height, bpp;
            string filePath = path+"/"+it->second->KdMap;
            unsigned char* tex = stbi_load(filePath.c_str(), &width, &height, &bpp, 3);
            vector<Vector> texture;
            string textureFile = it->first; 
            texture.resize(width * height);
            for (int i=0; i<width*height; i++){
                texture[i] = Vector(pow(tex[i*3] / 255., 2.2), pow(tex[i*3+1] / 255., 2.2), pow(tex[i*3+2] / 255., 2.2));
            }
            texture[0].print();
            it->second->texture = texture;
            it->second->width = width;
            it->second->height = height;
        }
    }
}