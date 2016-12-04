#ifndef _VECMAT_UTILS_
#define _VECMAT_UTILS_

extern int choldc1(double * a, double * p, int n);
extern int choldcsl(double * A, double * a, double * p, int n);
extern int cholsl(double * A, double * a, double * p, int n);
extern void zeros(double * a, int m, int n);
extern void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols);
extern void mulvec(double * a, double * x, double * y, int m, int n);
extern void transpose(double * a, double * at, int m, int n);
extern void accum(double * a, double * b, int m, int n);
extern void accum_scalar(double * a, double * b, int m, int n, double s1, double s2);
extern void add(double * a, double * b, double * c, int n);
extern void sub(double * a, double * b, double * c, int n);
extern void negate(double * a, int m, int n);
extern void mat_addeye(double * a, int n);

#endif 




