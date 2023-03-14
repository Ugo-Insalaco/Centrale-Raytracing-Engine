#ifndef BOXNODE_H
#define BOXNODE_H

#include "boundingBox.h"

class BoxNode {
    public:
        BoxNode(BoundingBox boundingBox, int startIndex, int endIndex);
        void printBoundaries();
        void print();
        BoundingBox boundingBox;
        BoxNode* leftChild {};
        BoxNode* rightChild {};
        int startIndex, endIndex;
};
#endif