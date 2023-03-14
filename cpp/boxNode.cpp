#include "../h/boxNode.h"
BoxNode::BoxNode(BoundingBox boundingBox, int startIndex, int endIndex) : boundingBox(boundingBox), startIndex(startIndex), endIndex(endIndex){};
void BoxNode::print(){
    cout << "BoxNode: startindex = " << startIndex << " endIndex = " << endIndex << endl;
    if(leftChild && rightChild){
        leftChild->print();
        rightChild->print();
    }
};
void BoxNode::printBoundaries(){
    cout << "BoxNode boundaries: startindex = " << startIndex << " endIndex = " << endIndex << endl;
}