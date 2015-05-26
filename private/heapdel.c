#include "mextools.h"
#include "heap.h"
#include "mex.h"

#define N_IN   prhs[0]
#define X_IN   prhs[1]
#define Y_IN   prhs[2]
#define N_OUT  plhs[0]
#define X_OUT  plhs[1]
#define Y_OUT  plhs[2]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int m, n;
    double *x, *y;

    if (nrhs != 3) {
        mexErrMsgTxt("Exactly 3 input arguments expected.");
        return;
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
        return;
    }
    
    m = mxGetM( X_IN );
    x = assertMatrix(X_IN, m, 1, "x");
    y = assertMatrix(Y_IN, m, 1, "y");
    n = (int)assertScalar(N_IN, "n"); 

    if (m < n) {
        mexErrMsgTxt("n is longer than size(x,1).");
        return;
    }

    if (nlhs >= 1) N_OUT = mxCreateDoubleScalar((double)(n-1));
    if (nlhs >= 2) X_OUT = mxCreateDoubleScalar(x[0]);
    if (nlhs >= 3) Y_OUT = mxCreateDoubleScalar(y[0]);

    heap_del_min(n,x,y);
}
