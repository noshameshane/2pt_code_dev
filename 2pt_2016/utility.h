/* Routines from Numerical Recipes (slightly modified). */
#ifndef UTILITY_H
#define UTILITY_H
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "constant.h"

using namespace std;

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define ludTINY 1.0e-20

int jacobi(double **a, int n, double d[], double **v);
int jacobi0(double **a, int n, double d[], double **v);
int eigsrt(double d[], double **v, int n);
int eigsrt0(double d[], double **v, int n);
int ROTATE(double **a,int i,int j,int k,int l,double *tau,double *s);
int solve_cubic_eq(double a,double b,double c,double *c3rts);
double scweighting(double upper);
double sqweighting(double upper);
double polylog(double n, double z);
double search2PT(double K);
double HSDF(double *pwr,int lenth,double fmin,int nmol,double fract_f,double fmf);
double TTDF(double *pwr,int lenth,double fmin);
int twoPTmf(double *pwr,double s0,int n,double pwrfreq,int nmol,double f2pt,double &fmf);
double cal_AB(double &A, double &B, double s0, double f2pt, double fmf);
double Sgmf(double v,double s0, double f2pt,int nmol,double fmf);
int twoPT(double *tmpg,double *tmps,double s0,double sv,double v,int nmol,double fract_f,double fmf);
int HSweighting(double *wep,double *wap,double *wcvp,double *wsehs,double *wsp,double *wspd,double y,double mass,double nmol,double tranT, double rotT,double volume,double *wsr,double *war,double *rT,double rs);
double dist(double [3],double[3]);
double angl(double [3],double[3],double[3]);
double dihe(double [3],double[3],double[3],double[3]);

int ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double *b);
int new_2d_matrix(double ** &M,int row,int col);
int delete_2d_matrix(double **&M,int row,int col);
int show_2d_matrix(double **&M,int row,int col);
int matrix_sum(double **A,double **B,int row,int col,double **C);
int matrix_sum(double **A,double **B,int row,int col,double **C,double fac);
int matrix_sub(double **A,double **B,int row,int col,double **C);
int matrix_product(double **A,int ra,int ca,double **B,int rb,int cb,double **C);
int matrix_product(double **A,int ra,int ca,double *B,int rb,double *C);
int matrix_BMBt(double **B,int rank,double **Mi,double **G); 
int matrix_inverse(double **A,int rank, double **invA, double &detA);
int matrix_scale(double **A,int ra,int ca,double scale);
int matrix_sqrt(double **A,int rank,double **B);

#endif
