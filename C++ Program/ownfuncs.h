#ifndef __ownfuncs_h__
#define __ownfuncs_h__

#include<vector>
#include<cmath>
#include<iostream>

int mod(int input, int ceil) {
    if (input > ceil) {
        return input - ceil;
    } else if (input < 0) {
        return input + ceil;
    } else {
        return input;
    }
}

double Vx(int i, int& L, int& D, std::vector<double>& sigma) {
    double Vx = 0;
    
    for (int d = 0; d < D; ++d) {
        Vx += sigma[mod( i - pow(L,d) , pow(L, D))];
        Vx += sigma[mod( i + pow(L,d) , pow(L, D))];
	}
    
    return Vx;
}

double VxHalf(int i, int& L, int& D, std::vector<double>& sigma) {
    double Vx = 0;
    
    for (int d = 0; d < D; ++d) {
        Vx += sigma[mod( i - pow(L,d) , pow(L, D))];
    }

    return Vx;
}

void hopping(int& L, int& D, std::vector<std::vector<int>>& hop) {
    int x, y, Lk;
    int xk, k, dxk ;

    /* go through all the points*/
    for (x=0; x < pow(L,D) ; x++){
	Lk = pow(L,D);
	y  = x;

	/* go through the components k*/
	for (k=D-1; k >= 0; k--){

	    Lk/=L;                        /* pow(L,k)      */
	    xk =y/Lk;                     /* kth component */
	    y  =y-xk*Lk;                  /* y<-y%Lk       */

	    /* forward */
	    if (xk<L-1) dxk = Lk;
	    else        dxk = Lk*(1-L);
	    hop[x][k] = x + dxk;

	    /* backward */
	    if (xk>0)   dxk = -Lk;
	    else        dxk = Lk*(L-1);
	    hop[x][k+D] = x + dxk;

	}
    }
} /* hopping */

#endif 