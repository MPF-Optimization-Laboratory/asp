/* $Id: heap.h 339 2008-07-15 18:33:57Z mpf $ */

#ifndef __HEAP_H__
#define __HEAP_H__

#define dswap(a,b) { double c; c = (a); (a) = (b); (b) = c; }

int  heap_del_min(int numElems, double x[], double y[]);
void heap_sift(int root, int lastChild, double x[], double y[]);
void heap_build(int n, double x[], double y[]);

#endif
