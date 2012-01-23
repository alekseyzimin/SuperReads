/* 
   Try to do the basic red-black tree stuff (no deletions,) assuming
   all the trees refer to the same array (that is, they share the same
   array for all their data, with their values interspersed.)
   More comments to follow.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <src2/redBlackTreesInsertOnly.h>

// INORDER-TREE-WALK (Only to print out in case of ints for now)
void inOrderTreeWalk (struct RBTreeStruct *T, int x, void (*func)(char *))
{
     void *vptr;

     if (x == TREE_NIL) return;
     inOrderTreeWalk (T, treeLeft (T, x), func);
     setTreeValPtr (vptr, T, x);
     
     func (vptr);
     inOrderTreeWalk (T, treeRight (T, x), func);
}

// RB-TREE-INSERT-IF-NOT-FOUND
void RBTreeInsertIfNotFound (struct RBTreeStruct *T, char *element)
{
     if (treeFindElement (T, element) == TREE_NIL)
	  RBTreeInsertElement (T, element);
}

// TREE-FIND-ELEMENT
int treeFindElement (struct RBTreeStruct *T, char *element)
{
     return (iterativeTreeSearch (T, T->root, element));
}

// ITERATIVE-TREE-SEARCH (x, k)
int iterativeTreeSearch (struct RBTreeStruct *T, int x, char *element)
{
     char *vptr;
     int ret;

     while (x != TREE_NIL) {
	  setTreeValPtr (vptr, T, x);
	  ret = ((int (*)()) ((T->dataArrayPtr)->compareToSearch)) (element, vptr);
	  if (ret == 0) break;
	  if (ret < 0)
	       x = treeLeft(T, x);
	  else
	       x = treeRight(T, x);
     }
     return (x);
}
	  
// TREE-MINIMUM (returns the index of the min value of the tree)
int treeMinimum (struct RBTreeStruct *T, int x)
{
     while (treeLeft(T, x) != TREE_NIL)
	  x = treeLeft(T, x);

     return x;
}

// TREE-MAXIMUM (returns the index of the max value of the tree)
int treeMaximum (struct RBTreeStruct *T, int x)
{
     while (treeRight(T, x) != TREE_NIL)
	  x = treeRight(T, x);

     return x;
}

// TREE-SUCCESSOR
// x is an index into the array holding the tree
// It returns an index into the array
int treeSuccessor (struct RBTreeStruct *T, int x)
{
     int y;

     if (treeRight(T, x) != TREE_NIL)
	  return (treeMinimum (T, treeRight(T, x)));
     y = treeParent(T, x);
     while ((y != TREE_NIL) && (x == treeRight(T, y))) {
	  x = y;
	  y = treeParent(T, y);
     }

     return (y);
}

// TREE-PREDECESSOR
// x is an index into the array holding the tree
// It returns an index into the array
int treePredecessor (struct RBTreeStruct *T, int x)
{
     int y;

     if (treeLeft(T, x) != TREE_NIL)
	  return (treeMaximum (T, treeLeft(T, x)));
     y = treeParent(T, x);
     while ((y != TREE_NIL) && (x == treeLeft(T, y))) {
	  x = y;
	  y = treeParent(T, y);
     }

     return (y);
}

// TREE-INSERT
// Note: node z must have already been created with data, but 
// without putting it in a tree
void treeInsert (struct RBTreeStruct *T, int z)
{
     int x, y;
     char *vptrx, *vptry, *vptrz;

     y = TREE_NIL;
     x = T->root;
     while (x != TREE_NIL) {
	  y = x;
	  setTreeValPtr (vptrz, T, z);
	  setTreeValPtr (vptrx, T, x);
	  if (((int (*)()) T->dataArrayPtr->compareToSort) (vptrz, vptrx) < 0)
	       x = treeLeft(T, x);
	  else
	       x = treeRight(T, x);
     }
     treeParent(T, z) = y;
     if (y == TREE_NIL) {
	  T->root = z;
	  return;
     }
     setTreeValPtr (vptry, T, y);
     setTreeValPtr (vptrz, T, z);
     if (((int (*)()) T->dataArrayPtr->compareToSort) (vptrz, vptry) < 0)
	  treeLeft(T, y) = z;
     else
	  treeRight(T, y) = z;

}

// TREE-LEFT-ROTATE
void treeLeftRotate (struct RBTreeStruct *T, int x)
{
     int y;

     y = treeRight (T, x);
     treeRight (T, x) = treeLeft (T, y);
     if (treeLeft (T, y) != TREE_NIL)
	  treeParent (T, treeLeft (T, y)) = x;
     treeParent (T, y) = treeParent (T, x);
     if (treeParent (T, x) == TREE_NIL) {
	  T->root = y;
	  goto skipTests1;
     }
     if (x == treeLeft (T, treeParent (T, x))) {
	  treeLeft (T, treeParent (T, x)) = y;
	  goto skipTests1;
     }
     treeRight (T, treeParent (T, x)) = y;

 skipTests1:
     treeLeft (T, y) = x;
     treeParent (T, x) = y;
}

// TREE-RIGHT-ROTATE
void treeRightRotate (struct RBTreeStruct *T, int x)
{
     int y;

     y = treeLeft (T, x);
     treeLeft (T, x) = treeRight (T, y);
     if (treeRight (T, y) != TREE_NIL)
	  treeParent (T, treeRight (T, y)) = x;
     treeParent (T, y) = treeParent (T, x);
     if (treeParent (T, x) == TREE_NIL) {
	  T->root = y;
	  goto skipTests2;
     }
     if (x == treeRight (T, treeParent (T, x))) {
	  treeRight (T, treeParent (T, x)) = y;
	  goto skipTests2;
     }
     treeLeft (T, treeParent (T, x)) = y;

 skipTests2:
     treeRight (T, y) = x;
     treeParent (T, x) = y;
}

// RB-TREE-INSERT-NODE
void RBTreeInsertNode (struct RBTreeStruct *T, int x)
{
     int y;

     treeInsert (T, x);
     treeColor (T, x) = RED;
     while ((x != T->root) && (treeColor (T, treeParent (T, x)) == RED)) {
	  if (treeParent (T, x) == treeLeft (T, treeParent (T, treeParent (T, x)))) {
	       y = treeRight (T, treeParent (T, treeParent (T, x)));
	       if (treeColor (T, y) == RED) {
		    treeColor (T, treeParent (T, x)) = BLACK;
		    treeColor (T, y) = BLACK;
		    treeColor (T, treeParent (T, treeParent (T, x))) = RED;
		    x = treeParent (T, treeParent (T, x));
	       }
	       else {
		    if (x == treeRight (T, treeParent (T, x))) {
			 x = treeParent (T, x);
			 treeLeftRotate (T, x);
		    }
		    treeColor (T, treeParent (T, x)) = BLACK;
		    treeColor (T, treeParent (T, treeParent (T, x))) = RED;
		    treeRightRotate (T, treeParent (T, treeParent (T, x)));
	       }
	  }
	  else {
	       y = treeLeft (T, treeParent (T, treeParent (T, x)));
	       if (treeColor (T, y) == RED) {
		    treeColor (T, treeParent (T, x)) = BLACK;
		    treeColor (T, y) = BLACK;
		    treeColor (T, treeParent (T, treeParent (T, x))) = RED;
		    x = treeParent (T, treeParent (T, x));
	       }
	       else {
		    if (x == treeLeft (T, treeParent (T, x))) {
			 x = treeParent (T, x);
			 treeRightRotate (T, x);
		    }
		    treeColor (T, treeParent (T, x)) = BLACK;
		    treeColor (T, treeParent (T, treeParent (T, x))) = RED;
		    treeLeftRotate (T, treeParent (T, treeParent (T, x)));
	       }
	  }
     }
     treeColor (T, T->root) = BLACK;
}

// RB-TREE-INSERT-ELEMENT
void RBTreeInsertElement (struct RBTreeStruct *T, char *element)
{
     struct dataArrayStruct *dasptr = T->dataArrayPtr;
     char *vptr;
     int numElemsHold, amountToAllocate;
     struct RBTreeNodeStruct *rbtptr;

//     fprintf (stderr, "Array size = %d\n", dasptr->arraySize); // DEBUGGING 1/20/12
     if (dasptr->arraySize == dasptr->numElementsAllocated) {
	  numElemsHold = dasptr->numElementsAllocated;
//	  printf ("numElemsHold = %d\n", numElemsHold);
	  (dasptr->numElementsAllocated) *= 2;
//	  printf ("dasptr->numElementsAllocated = %d\n", dasptr->numElementsAllocated);
//	  printf ("dasptr->elementSize = %d\n", dasptr->elementSize);
	  dasptr->nodeStructs = realloc (dasptr->nodeStructs, dasptr->numElementsAllocated * sizeof (struct RBTreeNodeStruct));
	  if (dasptr->nodeStructs == 0) {
	       fprintf (stderr, "realloc failed in RBTreeInsertElement(3). Bye!\n");
	       exit (3);
	  }
	  rbtptr = dasptr->nodeStructs + numElemsHold;
	  memset (rbtptr, 0, numElemsHold * sizeof (struct RBTreeNodeStruct));

	  amountToAllocate = dasptr->numElementsAllocated * dasptr->elementSize;
//          printf ("dasptr->numElementsAllocated = %d\n", dasptr->numElementsAllocated);
//          printf ("dasptr->elementSize = %d\n", dasptr->elementSize);
//	  printf ("amountToAllocate = %d\n", amountToAllocate);
//	  printf ("dasptr->data = %p\n", dasptr->data);
	  dasptr->data = realloc (dasptr->data, amountToAllocate);
	  if (dasptr->data == 0) {
	       fprintf (stderr, "realloc failed to allocate %d bytes in RBTreeInsertElement(4). numElements = %d, elementSize = %d. Bye!\n", amountToAllocate, dasptr->numElementsAllocated, dasptr->elementSize);
	       exit (4);
	  }
	  vptr = dasptr->data + (numElemsHold * dasptr->elementSize);
	  memset (vptr, 0, numElemsHold * dasptr->elementSize);
     }
     
     setTreeValPtr (vptr, T, dasptr->arraySize);
     memcpy (vptr, element, dasptr->elementSize);

     treeLeft (T, dasptr->arraySize) = TREE_NIL;
     treeRight (T, dasptr->arraySize) = TREE_NIL;

     dasptr->arraySize++;

     RBTreeInsertNode (T, (dasptr->arraySize) - 1);
}

