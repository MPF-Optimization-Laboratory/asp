#include <octave/oct.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include "heap.h"
#ifdef __cplusplus
}  /* end extern "C" */
#endif

DEFUN_DLD(heapify, args, nargout, "heapify input vectors")
{
    octave_value_list retval;

    if (args.length() != 2) {
        error("heapify: two input arguments expected.");
        return retval;
    }
    
    Matrix x = args(0).matrix_value();
    Matrix y = args(1).matrix_value();
    if (error_state) {
        error("heapify: problem with inputs");
        return octave_value_list();
    }
    int m = x.nelem();

    double *xv = x.fortran_vec();
    double *yv = y.fortran_vec();

    heap_build(m, xv, yv);

    retval(0) = x;
    retval(1) = y;
    return retval;

}
