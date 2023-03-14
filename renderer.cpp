#define _CRT_SECURE_NO_WARNINGS 1

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./h/stb_image_write.h"

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

#include <omp.h>
#include <cmath>
#include <iostream>
#include <stack>
#include <random>
#include <list>
#include <chrono>
#include <vector>
#include "./cpp/vector.cpp"
#include "./cpp/utils.cpp"
#include "./cpp/material.cpp"
#include "./cpp/ray.cpp"
#include "./h/geometry.h"
#include "./cpp/mesh.cpp"
#include "./cpp/sphere.cpp"
#include "./cpp/camera.cpp"
#include "./cpp/scene.cpp"

using namespace std;

int main() {
    #include "./scenes/mariokart.cpp"
    auto start = chrono::high_resolution_clock::now();
    

    mk_scene.render();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << duration.count() << endl;
    // double temp;
    // cin >> temp;
    return 0;
}