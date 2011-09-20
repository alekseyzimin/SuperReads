/* 
   Try to do the basic red-black tree stuff (no deletions,) assuming
   all the trees refer to the same array (that is, they share the same
   array for all their data, with their values interspersed.)
   More comments to follow.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TREE_DATA_START_SIZE 50000000

#define TREE_NIL -1
#define RED 'R'
#define BLACK 'B'

#define treeLeft(T, x) ((T->dataArrayPtr)->nodeStructs)[x].left
#define treeRight(T, x) ((T->dataArrayPtr)->nodeStructs)[x].right
#define treeParent(T, x) ((T->dataArrayPtr)->nodeStructs)[x].parent
#define treeColor(T, x) ((T->dataArrayPtr)->nodeStructs)[x].color

struct RBTreeStruct {
     int root;
     // The following is the ptr to the area holding the actual data + more
     struct dataArrayStruct *dataArrayPtr;
};

struct RBTreeNodeStruct {
     int parent;
     int left;
     int right;
     char color;
};

struct dataArrayStruct {
     int arraySize;
     int numElementsAllocated;
     int elementSize;
     struct RBTreeNodeStruct *nodeStructs; // The connection data
     void *data; // The actual data
     void *compareToSort;
     void *compareToSearch;
};

#define setTreeValPtr(ptr, T, ind) { ptr = ((T)->dataArrayPtr)->data; ptr += (ind * (((T)->dataArrayPtr)->elementSize)); }

#ifndef mallocOrDie
#define mallocOrDie(name, num, type) fprintf (stderr, "Allocating %lu bytes for %s.\n", (unsigned long) ((num) * sizeof ( type )), #name); \
name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }
#endif

#define initializeEmptyTrees(treeArrayName, numTrees, dataArrayInfo, dataType, sortCompare, searchCompare) { int i1234567; \
   dataArrayInfo.arraySize = 0; \
   dataArrayInfo.numElementsAllocated = TREE_DATA_START_SIZE; \
   dataArrayInfo.elementSize = sizeof (dataType); \
   mallocOrDie (dataArrayInfo.nodeStructs, TREE_DATA_START_SIZE, struct RBTreeNodeStruct); \
   mallocOrDie (dataArrayInfo.data, TREE_DATA_START_SIZE, dataType); \
   dataArrayInfo.compareToSort = sortCompare; \
   dataArrayInfo.compareToSearch = searchCompare; \
   mallocOrDie (treeArrayName, numTrees, struct RBTreeStruct); \
   for (i1234567=0; i1234567<numTrees; i1234567++) { \
      treeArrayName[i1234567].root = TREE_NIL; \
      treeArrayName[i1234567].dataArrayPtr = &dataArrayInfo; \
   } \
}
 
#define initializeEmptyTreesWithDataSize(treeArrayName, numTrees, dataArrayInfo, dataType, dataSize, sortCompare, searchCompare) { int i1234567; \
   dataArrayInfo.arraySize = 0; \
   dataArrayInfo.numElementsAllocated = dataSize; \
   dataArrayInfo.elementSize = sizeof (dataType); \
   mallocOrDie (dataArrayInfo.nodeStructs, dataSize, struct RBTreeNodeStruct); \
   mallocOrDie (dataArrayInfo.data, dataSize, dataType); \
   dataArrayInfo.compareToSort = sortCompare; \
   dataArrayInfo.compareToSearch = searchCompare; \
   mallocOrDie (treeArrayName, numTrees, struct RBTreeStruct); \
   for (i1234567=0; i1234567<numTrees; i1234567++) { \
      treeArrayName[i1234567].root = TREE_NIL; \
      treeArrayName[i1234567].dataArrayPtr = &dataArrayInfo; \
   } \
}
 
void inOrderTreeWalk (struct RBTreeStruct *T, int x, void (*)(void *));
void RBTreeInsertIfNotFound (struct RBTreeStruct *T, void *element);
int treeFindElement (struct RBTreeStruct *T, void *element);
int iterativeTreeSearch (struct RBTreeStruct *T, int x, void *element);
int treeMinimum (struct RBTreeStruct *T, int x);
int treeMaximum (struct RBTreeStruct *T, int x);
int treeSuccessor (struct RBTreeStruct *T, int x);
int treePredecessor (struct RBTreeStruct *T, int x);
void treeInsert (struct RBTreeStruct *T, int z);
void treeLeftRotate (struct RBTreeStruct *T, int x);
void treeRightRotate (struct RBTreeStruct *T, int x);
void RBTreeInsertNode (struct RBTreeStruct *T, int x);
void RBTreeInsertElement (struct RBTreeStruct *T, void *element);

