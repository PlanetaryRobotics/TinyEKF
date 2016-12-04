/*
 * TinyEKF: Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * This code is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ekf-> code.  If not, see <http:#www.gnu.org/licenses/>.
 */

/*
 * Change history
 * -------------
 * 11/28/2016  - Changes to make the state and state-msmt transitions configurable
 *             - Pull out matrix operations into a separate file (vecmat_utils)
 */

/* TODOs
 * -----
 * - Use Joseph's formula for computing Pk
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "vecmat_utils.h"

#include "tiny_ekf.h"

typedef struct 
{
    double * x;  /* state vector */

    double * A;  /* state transition function */
    double * C;  /* state measurement function */
    double * P;  /* prediction error covariance */
    double * Q;  /* process noise covariance */
    double * R;  /* measurement error covariance */

    double * G;  /* Kalman gain; a.k.a. K */

    double * F;  /* Jacobian of process model */
    double * H;  /* Jacobian of measurement model */

    double * Ht; /* transpose of measurement Jacobian */
    double * Ft; /* transpose of process Jacobian */
    double * Pp; /* P, post-prediction, pre-update */

    double * fx;  /* output of user defined f() state-transition function */
    double * hx;  /* output of user defined h() measurement function */

    /* temporary storage */
    double * tmp1;
    double * tmp2;
    double * tmp3;
    double * tmp4;
    double * tmp5; 

} ekf_t;

static void unpack(void * v, ekf_t * ekf, int n, int m)
{
    /* skip over n, m in data structure */
    char * cptr = (char *)v;
    cptr += 2*sizeof(int);

    double * dptr = (double *)cptr;
    ekf->x = dptr;
    dptr += n;
    ekf->A = dptr;
    dptr += n*n;
    ekf->C = dptr;
    dptr += m*n;
    ekf->P = dptr;
    dptr += n*n;
    ekf->Q = dptr;
    dptr += n*n;
    ekf->R = dptr;
    dptr += m*m;
    ekf->G = dptr;
    dptr += n*m;
    ekf->F = dptr;
    dptr += n*n;
    ekf->H = dptr;
    dptr += m*n;
    ekf->Ht = dptr;
    dptr += n*m;
    ekf->Ft = dptr;
    dptr += n*n;
    ekf->Pp = dptr;
    dptr += n*n;
    ekf->fx = dptr;
    dptr += n;
    ekf->hx = dptr;
    dptr += m;
    ekf->tmp1 = dptr;
    dptr += n*m;
    ekf->tmp2 = dptr;
    dptr += m*n;
    ekf->tmp3 = dptr;
    dptr += m*m;
    ekf->tmp4 = dptr;
    dptr += m*m;
    ekf->tmp5 = dptr;
  }

void ekf_init(void * v, int n, int m)
{
    /* retrieve n, m and set them in incoming data structure */
    int * ptr = (int *)v;
    *ptr = n;
    ptr++;
    *ptr = m;

    /* unpack rest of incoming structure for initlization */
    ekf_t ekf;
    unpack(v, &ekf, n, m);

    /* zero-out matrices */
    zeros(ekf.A, n, n);
    zeros(ekf.C, m, n);
    zeros(ekf.P, n, n);
    zeros(ekf.Q, n, n);
    zeros(ekf.R, m, m);
    zeros(ekf.G, n, m);
    zeros(ekf.F, n, n);
    zeros(ekf.H, m, n);
}

void ekf_predict(void * v, double delta)
{
    /* unpack incoming structure */

    int * ptr = (int *)v;
    int n = *ptr;
    ptr++;
    int m = *ptr;

    ekf_t ekf;
    unpack(v, &ekf, n, m);

    /* Predict step */
    mulmat(ekf.A, ekf.x, ekf.fx, n, n, 1);		       

    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
    mulmat(ekf.F, ekf.P, ekf.tmp1, n, n, n);
    transpose(ekf.F, ekf.Ft, n, n);
    mulmat(ekf.tmp1, ekf.Ft, ekf.Pp, n, n, n);
    accum_scalar(ekf.Pp, ekf.Q, n, n, 1.0, delta);
}

int ekf_step(void * v, double * z, int idx)
{        
    /* unpack incoming structure */
    int * ptr = (int *)v;
    int n = *ptr;
    ptr++;
    int m = *ptr;

    ekf_t ekf;
    unpack(v, &ekf, n, m);
		   
    /* hx = C*fx */
    mulvec(ekf.C, ekf.fx, ekf.hx, m, n);		           
    
    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
    transpose(ekf.H, ekf.Ht, m, n);
    mulmat(ekf.Pp, ekf.Ht, ekf.tmp1, n, n, m);
    mulmat(ekf.H, ekf.Pp, ekf.tmp2, m, n, n);
    mulmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, n, m);
    accum(ekf.tmp3, ekf.R, m, m);
    if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m)) return 1;
    mulmat(ekf.tmp1, ekf.tmp4, ekf.G, n, m, m);

    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
    sub(z, ekf.hx, ekf.tmp5, m);
    mulvec(ekf.G, ekf.tmp5, ekf.tmp2, n, m);
    add(ekf.fx, ekf.tmp2, ekf.x, n);

    /* P_k = (I - G_k H_k) P_k */ 
    mulmat(ekf.G, ekf.H, ekf.tmp1, n, m, n);
    negate(ekf.tmp1, n, n);
    mat_addeye(ekf.tmp1, n);
    mulmat(ekf.tmp1, ekf.Pp, ekf.P, n, n, n);

    /* success */
    return 0;
}
