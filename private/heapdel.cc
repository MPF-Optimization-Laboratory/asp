#include <octave/oct.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include "heap.h"
#ifdef __cplusplus
}  /* end extern "C" */
#endif

DEFUN_DLD(heapdel, args, nargout, "Delete smallest element from heap")
{
    octave_value_list retval;
    int nrhs = args.length();

    if (nrhs != 3 || nrhs > 3 ) {
        error("heapdel: expecte 3 input and 5 outputs.");
        return retval;
    }

    Matrix x = args(1).matrix_value();
    Matrix y = args(2).matrix_value();
    int n = args(0).int_value();
    int m = x.nelem();
    double *xv = x.fortran_vec();
    double *yv = y.fortran_vec();

    if (m < n) {
        error("heapdel: n is longer than size(x,1).");
        return retval;
    }

    retval(0) = octave_value(n-1);
    retval(1) = octave_value(xv[0]);
    retval(2) = octave_value(yv[0]);
    retval(3) = x;
    retval(4) = y;

    heap_del_min(n,xv,yv);

    return retval;
}
