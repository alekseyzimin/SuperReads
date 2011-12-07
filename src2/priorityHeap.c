/*
  Try to do the basic heap stuff
  The only stuff not to be ported to another application is guts of main
  and the declaration and definition of intCompare.
*/
#define PRIORITY_HEAP_C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <src2/priorityHeap.h>

int heapIsEmpty (struct heapStruct *A)
{
     if (A->heapSize == 0) return 1;
     else                  return 0;
}

void heapInsert (struct heapStruct *A, void *element)
{
     void *vptr1, *vptr2;
     int i;

     if (A->heapSize == A->numElementsAllocated) {
	  (A->numElementsAllocated) *= 2;
	  A->data = realloc (A->data, (A->numElementsAllocated) * (A->elementSize));
	  if (A->data == NULL) {
	       fprintf (stderr, "Ran out of data for realloc. Bye!\n");
	       exit (2);
	  }
     }
     (A->heapSize)++;
     i = A->heapSize - 1;
     while (i > 0) {
	  setHeapValPtr (vptr1, A, HEAP_PARENT(i));
	  if (((int (*)()) A->compare) (vptr1, element) < 0) {
	       setHeapValPtr (vptr2, A, i);
	       memcpy (vptr2, vptr1, A->elementSize);
	       i = HEAP_PARENT (i);
	  }
	  else
	       break;
     }
     setHeapValPtr(vptr1, A, i);
     memcpy (vptr1, element, A->elementSize);
}

// Note that the space for element must exist (i.e. call as &element)
void heapExtractRoot (struct heapStruct *A, void *element)
{
     void *vptr;

     if (A->heapSize < 1) {
	  fprintf (stderr, "No elements left to extract. Bye!\n");
	  exit (1);
     }
     memcpy (element, A->data, A->elementSize);
     (A->heapSize)--;
     if (A->heapSize > 0) {
	  setHeapValPtr (vptr, A, A->heapSize);
	  memcpy (A->data, vptr, A->elementSize);
     }
     heapify (A, 0);
}

void heapify (struct heapStruct *A, int i)
{
     int l, r, largest;
     unsigned char *vptr1, *vptr2;

     l = HEAP_LEFT (i);
     r = HEAP_RIGHT (i);
     largest = -1;
     if (l < A->heapSize) {
	  setHeapValPtr (vptr1, A, l);
	  setHeapValPtr (vptr2, A, i);
	  if (((int (*)()) A->compare) (vptr1, vptr2) > 0)
	       largest = l;
     }
     if (largest < 0) 
	  largest = i;
     if (r < A->heapSize) {
	  setHeapValPtr (vptr1, A, r);
	  setHeapValPtr (vptr2, A, largest);
	  if (((int (*)()) A->compare) (vptr1, vptr2) > 0)
	       largest = r;
     }
     if (largest != i) {
	  // exchange A[i] and A[largest]
	  setHeapValPtr (vptr1, A, i);
	  setHeapValPtr (vptr2, A, largest);
	  memcpy (tempHeapPtrs[A->tempSpaceNumber], vptr1, A->elementSize);
	  memcpy (vptr1, vptr2, A->elementSize);
	  memcpy (vptr2, tempHeapPtrs[A->tempSpaceNumber], A->elementSize);
	  heapify (A, largest);
     }
}

