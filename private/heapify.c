#include "heap.h"
#include "mextools.h"
#include "mex.h"

#define X_IN prhs[0]
#define Y_IN prhs[1]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int m;
    double *x, *y;

    if (nrhs != 2) {
        mexErrMsgTxt("Two input arguments expected.");
        return;
    }

    m = mxGetM( X_IN );
    x = assertMatrix(X_IN, m, 1, "x");
    y = assertMatrix(Y_IN, m, 1, "y");

    heap_build(m,x,y);
}
