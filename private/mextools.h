/* $Id: mextools.h 339 2008-07-15 18:33:57Z mpf $ */
/* Matlab helper routines. */

#ifndef __MEXTOOLS_H__
#define __MEXTOOLS_H__

#include "mex.h"

typedef struct {  /* Matlab full matrix. */
  int m, n;
  double *val;
} matlab_Mat;

typedef struct {  /* Compressed column sparse format. */
  int    n, m, nnz;
  int    *ind, *col;
  double *val;
} matlab_spMat;

static double
assertScalar( const mxArray *mex, char who[8] )
{
   char msgbuf[60];
   int  m = mxGetM( mex );
   int  n = mxGetN( mex );
   if ( m != 1 || n != 1 || mxIsSparse( mex ) ) {
      sprintf( msgbuf,
	       "Expected '%.8s' to be a scalar.  Instead it's %dx%d.\n",
	       who,m,n);
      mexErrMsgTxt( msgbuf );
   }
   else {
      return mxGetScalar( mex );
   }
}

static char *
assertString( const mxArray *mex, char who[8], int add )
{
  char msgbuf[60];
  int  m = mxGetM( mex );
  int  n = mxGetN( mex );
  if ( m != 1 || mxIsSparse( mex ) ) {
    sprintf( msgbuf,
	     "Expected '%.8s' to be a single string.  Instead it's %dx%d.\n",
	     who,m,n);
    mexErrMsgTxt( msgbuf );
  }
  else {
    int  buflen  = mxGetN( mex ) + 1;
    char *string = mxCalloc( buflen + add, sizeof(char) );
    int err      = mxGetString( mex, string, buflen );
    if (err) {
      sprintf( msgbuf, "Not enough space.  '%.8s' was trucated.\n", who );
      mexWarnMsgTxt( msgbuf );
    }
    else {
      return string;
    }
  }
}

static double *
assertSpMatrix( const mxArray *mex,
		int m, int n, char who[8] )
{
   char msgbuf[60];
   int  m_actual = mxGetM( mex );
   int  n_actual = mxGetN( mex );
   if ( m != m_actual || n != n_actual || !mxIsSparse( mex ) ) {
      sprintf( msgbuf,
       "Expected '%.8s' to be a sparse %dx%d matrix.  Instead it's %dx%d.\n",
       who, m,n, m_actual,n_actual);
      mexErrMsgTxt( msgbuf );
   }
   else {
      return mxGetPr( mex );
   }
}

static double *
assertMatrix( const mxArray *mex,
	      int m, int n, char who[8] )
{
   char msgbuf[60];
   int  m_actual = mxGetM( mex );
   int  n_actual = mxGetN( mex );
   if ( m != m_actual || n != n_actual || mxIsSparse( mex ) ) {
      sprintf( msgbuf,
       "Expected '%.8s' to be a dense %dx%d matrix.  Instead it's %dx%d.\n",
       who, m,n, m_actual,n_actual);
      mexErrMsgTxt( msgbuf );
   }
   else {
      return mxGetPr( mex );
   }
}

static matlab_spMat
create_matlab_spmat( mxArray **plhs, int m, int n, int nnz )
{
  matlab_spMat A;
  A.n   = n;
  A.m   = m;
  A.nnz = nnz;
  *plhs = mxCreateSparse( m, n, nnz, mxREAL );
  A.val =       mxGetPr( *plhs );
  A.ind = (int*)mxGetIr( *plhs );
  A.col = (int*)mxGetJc( *plhs );
  return A;
}

#endif
