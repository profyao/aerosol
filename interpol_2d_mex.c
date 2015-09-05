#include "mex.h"
#include<stdio.h>
#include<stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  /*lookup point to 13*n array*/

  double tau_grid[13]={0,0.05,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3,4,6};
  double *tau, *lookup, *interpo, yy1, yy2, xx1=0, xx2=6;
  int    tau_m, tau_n, lookup_m, lookup_n, ini=0, inj, pos;

  if (nrhs != 2 || nlhs != 1)
    mexErrMsgTxt("Usage: inter = interpol_mex(tau, lookups)");

  tau = mxGetPr(prhs[0]);
  tau_m = mxGetM(prhs[0]);
  tau_n = mxGetN(prhs[0]);
  if(tau_m != 1 || tau_n != 1)
    mexErrMsgTxt("Usage: tau needs to be a scalar");

  lookup = mxGetPr(prhs[1]);
  lookup_m = mxGetM(prhs[1]);
  lookup_n = mxGetN(prhs[1]);
  /*if(lookup_m != 13 || lookup_n != 9)*/
  /*  mexErrMsgTxt("Usage: lookup needs to be a 13x9.");*/

  plhs[0] = mxCreateDoubleMatrix(lookup_n, 1, mxREAL);
  interpo = mxGetPr(plhs[0]);

  if ((*tau)<0){
    interpo[0] = -1; /*error message*/
	mexPrintf("tau == %f \n", *tau);}
  else if ((*tau)>6)
    for (ini = 0; ini<lookup_n; ini++)
      interpo[ini] = *(lookup+12+ini*13);
  else {
    xx1 = tau_grid[0];/*left bracket*/
    inj = 1;
    while (!((xx1<(*tau)) && ((*tau)<=tau_grid[inj]))) 
      xx1 = tau_grid[inj++];
    xx2 = tau_grid[inj];/*right bracket*/
    for (ini = 0; ini < lookup_n; ini++){
      pos = ini*13;
      yy1 = *(lookup+pos+inj-1);
      yy2 = *(lookup+pos+inj);
      interpo[ini] = yy1 + ((*tau)-xx1)*(yy2-yy1)/(xx2-xx1);
    }
  }
}
