#pragma once
#include "math.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>

#include <vector>
using std::vector;

#define DEPS 2.22044604925031e-016
// The one and only application object
#define R_MAX_Iter 0x01
#define R_MIN_TolFun 0x02
#define R_MIN_TolX 0x03
#define R_INC_SSE 0x04
#define R_ERROR 0x05
#define MAX_DOUBLE (10e38)

#define MX_ITER_LIM 100


#define PI 3.14159265358979
#define EPS 2e-52
#define TUNE 2.986

struct gsl_Data 
{
	int n;
	double *y;
	double *x;
};

struct guass_data
{
	size_t fn;
	size_t n;
	double *x;
	double *y;
	double *weight;
};

typedef struct gsl_Data FIT_DATA;

typedef double(*F)(double x,double *beta);
typedef double(*FNew)(double x,double *beta,int num);

struct gsl_nlfit_options
{
	int MaxIter;
	double TolFun;
	double TolX;
	double DerivStep;
	int Disp;
	double sse;
	double r;
};

typedef gsl_nlfit_options OP;
	///////////////////////////////////////////////////////////////////////////////////
double Gauss_F(double x, double *beta);
double gsl_gauss(double x, void *params);
void gsl_diff(double *x,double *diff,int n);
double gsl_pow(double x,int j);
bool sgolay(double *x,double *y,int n,int f,int k);
int gsl_fit_max(double *y,int n);
double mgsl_fit_min(double *y,int n);
bool Initial_c(FIT_DATA *d,double c[4]);
double gsl_norm(double *v,int n);
double gsl_sum(double *v,int n);
int gsl_divide(gsl_matrix *A,double *b,double *x);
double gsl_update_beta(gsl_matrix *J,double lambda,double *rplus,double *beta);
int gsl_nlfit(double *x,double *y,size_t n,double *beta,size_t p,F f,OP *options);
double Gauss_Fit(FIT_DATA *d,double beta[4]);
double GetArea(double beta[3],double a,double b);
gsl_matrix *gsl_matrix_alloc_s(size_t row,size_t col);
bool gsl_GetPShape(double *x,double *y,int n,double *PH,double *PW,double *PA);
double gsl_CShape(double beta[4],double x);
double gsl_ChroFit(double *x,double *y,int n,double beta[5]);
double gsl_CShapeNew(double x,double beta[4]);

double gsl_GetPA(double *x,double *y,int n);
double gsl_GetPCA(double *li,double *hi,int n);
bool gsl_PSmooth(double *x,double *y,int n,int nnew);
double EM_Gauss_Fit(FIT_DATA *d,double beta[4],double L,double H);

int gsl_vector_output(gsl_vector *X);
int gsl_matrix_output(gsl_matrix *X);
int gsl_matrix_inverse(gsl_matrix *X);
double gsl_vector_sum(gsl_vector *X);
int gsl_matrix_mul(gsl_matrix *X,gsl_matrix *Y);
int gsl_matrix_div(gsl_matrix *X,gsl_matrix *Y);
bool gsl_adfactor_get(gsl_matrix *X,gsl_vector *h);
int gsl_sort_double(double *data,int n);
double gsl_stats_median(double *data,int n);
double gsl_madsigma(gsl_vector *Res,int p);
int gsl_update_w(gsl_vector *w,gsl_vector *Res,gsl_vector *h,int type,int p);
int gsl_res_get(gsl_matrix *X,gsl_vector *y,gsl_vector *c,gsl_vector *Res);
int gsl_res_get(gsl_matrix *X,gsl_vector *y,gsl_vector *c,double *Res);
int gsl_robustfit(gsl_matrix *X,gsl_vector *y, gsl_vector *c);
//double gsl_stats_median(double *data,int n);

double w_std(double *Res,char *ResW,int n,double m);
double w_mean(double *Res,char *ResW,int n);
double w_std(gsl_vector *Res,char *ResW,int n,double m);
double w_mean(gsl_vector *Res,char *ResW,int n);

double gsl_w_mean(gsl_vector *x,gsl_vector *w);
double gsl_mean(gsl_vector *x);
double gsl_w_std(gsl_vector *x,gsl_vector *w,double m);
//cluster fitting
double Gauss_F(double x, double *beta,int num);
double Gauss_F(double x, vector<double> &beta);
bool Initial_c(vector<double> &x,vector<double> &y,vector<double> &beta);
double Gauss_Fit(vector<double> &x,vector<double> &y,vector<double> &beta);

int gsl_nlfitNew(double *x,double *y,size_t n,double *beta,size_t p,FNew f,OP *options);
void Re_initial(vector<double> &x,vector<double> &y,vector<double> &beta);
//try orbitrap
bool Initial_cNew(vector<double> &x,vector<double> &y,vector<double> &beta);

int guass_f(const gsl_vector * par, void *data,gsl_vector * f);
int guass_df(const gsl_vector * par, void *data,gsl_matrix * J);
int guass_fdf (const gsl_vector * x, void *data,gsl_vector * f, gsl_matrix * J);
int gsl_nlfit_og(double *x,double *y,const size_t n,double *beta,size_t fn,gsl_nlfit_options *op);
double Gauss_F_new(double x, vector<double> &beta);
double Gauss_F_new(double x,double *beta,int num);
double Gauss_Fit_new(vector<double> &x,vector<double> &y,vector<double> &beta);

int weight_centriod(vector<double> &x,vector<double> &y);

double gsl_stats_correlation (const double data1[], const size_t stride1,const double data2[], const size_t stride2, const size_t n);

double mygsl_stats_mean(vector<double> &ppme);
double mygsl_stats_stdm(vector<double> &ppme,double m);
void mygsl_stats_stdm(vector<double> &ppme,double m[2]);
void mygsl_stats_stdm(vector<double> &ppme,double *m,double *sd);

double gsl_two_step_reg(double *x,double *y,int p,int n,double par[]);
int gsl_count(double *x,double lmin,double lmax,int n);
double gsl_normpdf(double x,double mu,double sigma);
double gsl_RSquare(gsl_vector *y,gsl_vector *Res,gsl_vector *w);

int gsl_robustfit(gsl_matrix *X,gsl_vector *y, gsl_vector *c,double stats[3]);