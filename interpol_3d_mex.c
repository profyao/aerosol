#include "mex.h"
#include<stdio.h>
#include<stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  /*lookup point to 13*9*k array*/

	const mwSize *lookup_dn;
	mwSize lookup_d;
	double tau_grid[13]={0,0.05,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3,4,6};
	double *tau, *lookup, *interpo, yy1, yy2, xx1=0, xx2=6, zz1, zz2;
	int    tau_m, tau_n, ini=0, inj, ink=0, pos, jj;


  if (nrhs != 2 || nlhs != 1)
    mexErrMsgTxt("Usage: inter = interpol_mex(tau, lookups)");

  tau = mxGetPr(prhs[0]);
  tau_m = mxGetM(prhs[0]);
  tau_n = mxGetN(prhs[0]);
  if(tau_m != 1 || tau_n != 1)
    mexErrMsgTxt("Usage: tau needs to be a scalar");

  lookup = mxGetPr(prhs[1]);
  lookup_d = mxGetNumberOfDimensions(prhs[1]);
  if(lookup_d != 3)
    mexErrMsgTxt("Usage: look-up table needs to be 3 dimensional");
  lookup_dn = mxGetDimensions(prhs[1]);
  /*if(lookup_dn[0] != 13 || lookup_dn[1] != 9 || lookup_dn[2] != 4)*/
  /*  mexErrMsgTxt("Usage: look-up table needs to be a 13x9x4 matrix")*/

  plhs[0] = mxCreateDoubleMatrix(9, lookup_dn[2], mxREAL);
  interpo = mxGetPr(plhs[0]);

  if ((*tau)<0){
    interpo[0] = -1; /*error message*/
	mexPrintf("Error: tau == %f \n", *tau);}
  else if ((*tau)>6)
	for (ink = 0; ink<lookup_dn[2]; ink++)
    	for (ini = 0; ini<9; ini++)
      		interpo[ini+ink*9] = *(lookup+12+ini*13+ink*117);/*need to check!!*/
  else {
    inj = 0;
    while ((*tau)>tau_grid[++inj]);
	inj--;
    xx1 = tau_grid[inj];/*left bracket*/
    xx2 = tau_grid[inj+1];/*right bracket*/
    switch (inj){
    	case 0:
    	case 1:
    	case 2:
    	case 3:
		case 9:
    	case 10:
    	case 11: /*linear interpolation*/
      		for (ink = 0; ink < lookup_dn[2]; ink++)
    			for (ini = 0; ini < 9; ini++){
      				pos = ini*13+ink*117;
      				yy1 = *(lookup+pos+inj);
      				yy2 = *(lookup+pos+inj+1);
      				interpo[ini+ink*9] = yy1 + ((*tau)-xx1)*(yy2-yy1)/(xx2-xx1);/*need to check!!!*/
    			}
			break;
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:/*quadratic interpolation*/
			for (ink = 0; ink < lookup_dn[2]; ink++)
    			for (ini = 0; ini < 9; ini++){
	      			pos = ini*13+ink*117;
					/*calculate all z_j values, j=4,...,inj+1 */
					zz1 = (*(lookup+pos+4)-*(lookup+pos+3))/(0.4-0.2);
					zz2 = -zz1 + 2*(*(lookup+pos+5)-*(lookup+pos+4))/(0.6-0.4);
					for (jj = 5; jj < inj+1; jj++){
						zz1 = zz2;
						zz2 = -zz1 + 2*(*(lookup+pos+jj+1)-*(lookup+pos+jj))/(tau_grid[jj+1]-tau_grid[jj]);
					}			
      				yy1 = *(lookup+pos+inj);
      				interpo[ini+ink*9] = yy1 + zz1*((*tau)-xx1) + ((*tau)-xx1)*((*tau)-xx1)*(zz2-zz1)/(2*(xx2-xx1));/*need to check!!!*/
    			}
			break;
    }
  }
}
