/*
  Try to do the basic heap stuff
  The only stuff not to be ported to another application is guts of main
  and the declaration and definition of intCompare.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HEAP_START_SIZE 50000000

#define HEAP_PARENT(i) (((i)-1)/2)
#define HEAP_LEFT(i) (2*(i)+1)
#define HEAP_RIGHT(i) (2*(i)+2)

struct heapStruct {
     int heapSize;
     int numElementsAllocated;
     void *compare;
     int elementSize;
     int tempSpaceNumber;
     void *data;
};

// The following holds temp space, one for each heap
#ifndef PRIORITY_HEAP_C
extern void *tempHeapPtrs[5000];
extern int numHeapsAllocated;
#else
void *tempHeapPtrs[5000];
int numHeapsAllocated=0;
#endif

#ifndef mallocOrDie
#define mallocOrDie(name, num, type) fprintf (stderr, "Allocating %lu bytes for %s.\n", (unsigned long) ((num) * sizeof ( type )), #name); \
name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }
#endif

#define initializeEmptyHeap(heapName, dataType, compareFunc) { mallocOrDie(heapName.data, HEAP_START_SIZE, dataType); \
heapName.heapSize = 0; \
heapName.numElementsAllocated = HEAP_START_SIZE; \
mallocOrDie (tempHeapPtrs[numHeapsAllocated], 1, dataType); \
heapName.tempSpaceNumber = numHeapsAllocated; \
numHeapsAllocated++; \
heapName.compare = (void *) compareFunc; \
heapName.elementSize = sizeof (dataType); }

#define setHeapValPtr(ptr, A, l) { ptr = (A)->data; ptr += (l * (A)->elementSize); }
void heapInsert (struct heapStruct *A, void *element);
void heapExtractRoot (struct heapStruct *A, void *element);
void heapify (struct heapStruct *A, int i);
int heapIsEmpty (struct heapStruct *A);

