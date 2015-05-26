#include "heap.h"

/*!
  
  \brief Discard the smallest element and contract the heap.

  On entry, the numElems of the heap are stored in x[0],...,x[numElems-1],
  and the smallest element is x[0].  The following operations are performed:
    -# Swap the first and last elements of the heap
    -# Shorten the length of the heap by one.
    -# Restore the heap property to the contracted heap.
       This effectively makes x[0] the next smallest element
       in the list.  

  \param[in]     numElems   The number of elements in the current heap.
  \param[in,out] x          The array to be heapified.
  \param[in,out] y          This array is correspondingly permuted.

  \return  The number of elements in the heap after it has been contracted.

*/
int
heap_del_min( int numElems, double x[], double y[] )
{
    int lastChild = numElems - 1;

    /* Swap the smallest element with the lastChild. */
    dswap( x[0], x[lastChild] );
    dswap( y[0], y[lastChild] );

    /* Contract the heap size, thereby discarding the smallest element. */
    lastChild--;
    
    /* Restore the heap property of the contracted heap. */
    heap_sift( 0, lastChild, x, y );

    return numElems - 1;
}

/*!

  \brief Perform the "sift" operation for the heap-sort algorithm.

  A heap is a collection of items arranged in a binary tree.  Each
  child node is less than or equal to its parent.  If x[k] is the
  parent, than its children are x[2k+1] and x[2k+2].

  This routine promotes ("sifts up") children that are smaller than
  their parents.  Thus, this is a "reverse" heap, where the smallest
  element of the heap is the root node.

  \param[in]     root       The root index from which to start sifting.
  \param[in]     lastChild  The last child (largest node index) in the sift operation.
  \param[in,out] x          The array to be sifted.
  \param[in,out] y          This array is correspondingly permuted.

*/
void
heap_sift( int root, int lastChild, double x[], double y[] )
{
    int child;

    for (; (child = (root * 2) + 1) <= lastChild; root = child) {

	if (child < lastChild)
	    if ( x[child] > x[child+1] )
		child++;
	
	if ( x[child] >= x[root] )
	    break;

	dswap( x[root], x[child] );
	dswap( y[root], y[child] );
    }
}

/*!
  
  \brief  Build a heap by adding one element at a time.
  
  \param[in]      n   The length of x and ix.
  \param[in,out]  x   The array to be heapified.
  \param[in,out]  y   This array is correspondingly permuted.

*/
void
heap_build( int n, double x[], double y[] )
{    
    int i;

    for (i = n/2; i >= 0; i--)
	heap_sift( i, n-1, x, y );
}
