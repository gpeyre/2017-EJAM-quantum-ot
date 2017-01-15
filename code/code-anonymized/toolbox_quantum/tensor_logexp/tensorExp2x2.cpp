#include <iostream>
// #include <math.h>
#include <cmath>
#include <map>
#include <cstring>
#include <vector>
#include "mex.h"
using namespace std;

extern void _main();

// vertices, faces, vectors, bases, radius, nRadialSamples
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray *prhs[]) {
	if (nrhs != 1)
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "TENSOREXP2X2 requires one input.");
    else if (nlhs != 1) // vertices, faces, indicator
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "TENSOREXP2X2 produces one output.");

    size_t nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd < 2) mexErrMsgTxt("TENSOREXP2X2 requires at least 2 array dimensions.");

    const mwSize *dim = mxGetDimensions(prhs[0]);
    if (dim[0] != 2 || dim[1] != 2)
        mexErrMsgTxt("TENSOREXP2X2 only works for 2x2 matrices.");

    double *data = mxGetPr(prhs[0]);

    size_t total = mxGetNumberOfElements(prhs[0]);
    if (total % 4 != 0)
        mexErrMsgTxt("TENSOREXP2X2 needs a total size divisible by 4.");

    size_t nTensors = total/4;

    double *result = (double*)mxCalloc(4*nTensors,sizeof(double));
    plhs[0] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
    mxSetData(plhs[0],result);
    mxSetDimensions(plhs[0],dim,nd);

    //#pragma omp parallel for
    for (int i = 0; i < nTensors; i++) {
        const double &a = data[4*i],
                     &b = data[4*i+3],
                     &c = data[4*i+2],
                     &d = data[4*i+1];

        if (abs(c-d) > 1e-9) mexErrMsgTxt("TENSORLOG2X2 only works for symmetric matrices.");

        double D = sqrt( (a-b)*(a-b) + 4*c*c ) / 2,
               T = (a+b) / 2,
               sign = ((a-b)>=0)? 1. : -1.,
               u = exp(T+sign*D),
               v = exp(T-sign*D);

        double theta = atan(2*c / (a-b))/2; // not atan2?

        if (isnan(theta)) theta = 0;

        double C = cos(theta), S = -sin(theta);
        double C2 = C*C;
        double S2 = S*S;
        double CS = C*S;

        result[4*i] = u*C2 + v*S2;
        result[4*i+1] = (-u+v)*CS;
        result[4*i+2] = result[4*i+1];
        result[4*i+3] = u*S2+v*C2;
    }
}
