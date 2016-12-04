/*
 * tiny_ekf_struct.h: common data structure for TinyEKF
 *
 * You should #include this file after using #define for N (states) and M
*  (observations)
 *
 * Copyright (C) 2016 Simon D. Levy
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
 * along with this code.  If not, see <http:#www.gnu.org/licenses/>.
 */

typedef struct 
{
    int n;          /* number of state values */
    int m;          /* number of observables */

    double x[N_STATE];    /* state vector */

    double A[N_STATE][N_STATE]; /* State-transition function */
    double C[N_MSMT][N_STATE];   /* State-measurement function */
    double P[N_STATE][N_STATE];  /* prediction error covariance */
    double Q[N_STATE][N_STATE];  /* process noise covariance */
    double R[N_MSMT][N_MSMT];    /* measurement error covariance */

    double G[N_STATE][N_MSMT];   /* Kalman gain; a.k.a. K */

    double F[N_STATE][N_STATE];  /* Jacobian of process model */
    double H[N_MSMT][N_STATE];   /* Jacobian of measurement model */

    double Ht[N_STATE][N_MSMT];  /* transpose of measurement Jacobian */
    double Ft[N_STATE][N_STATE]; /* transpose of process Jacobian */
    double Pp[N_STATE][N_STATE]; /* P, post-prediction, pre-update */

    double fx[N_STATE];   /* output of user defined f() state-transition function */
    double hx[N_MSMT];    /* output of user defined h() measurement function */

    /* temporary storage */
    double tmp1[N_STATE][N_STATE];
    double tmp2[N_MSMT][N_STATE];
    double tmp3[N_MSMT][N_MSMT];
    double tmp4[N_MSMT][N_MSMT];
    double tmp5[N_MSMT]; 

} ekf_t;        
