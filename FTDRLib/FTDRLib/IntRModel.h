#pragma once
#include <gsl/gsl_vector.h>
#include "gsl/gsl_multimin.h"
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>

#include <vector>
using std::vector;
#define TUNE 2.986
#define SQRTRPI 0.398942280401433
#define LOG2PIH -0.918938533204673
#define W_TUNE 4.0

class FittinfPar
{
public:
	vector<double> epsi;
	vector<double> yi;	
	vector<int> weight;
	double a;
	double b;
	double theta;
	FittinfPar();
	~FittinfPar();
	FittinfPar(const FittinfPar & r);
	FittinfPar &operator =(const FittinfPar &r);
	void Initial();
	double Error_mean();
	void UpdateW();
	double ParError(double aold,double bold);
};

double NormalPDF(double x,double mu,double sigma);
double LogNormalPDF(double x,double mu,double sigma);
double Error_mean(FittinfPar *Par);
double Error_std(FittinfPar *Par);
double IntR_Res_f(const gsl_vector *v, void *params);
void IntR_Res_df(const gsl_vector *v, void *params,gsl_vector *df);
void IntR_Res_fdf(const gsl_vector *v, void *params,double *f,gsl_vector *df);
double Fitting(FittinfPar *Par);
double Fitting_W(FittinfPar *Par);
bool IntR_ModelBuilding(double *ei,double *Ii,int n,double rbPar[3]);
bool IntR_ModelBuilding(vector<double> &ei,vector<double> &Ii,double rbPar[3]);
bool OutPutModelData(char *fname,vector<double> &ei,vector<double> &Ii);
bool OutPutModelData(char *fname,double *ei,double *Ii,int n);

bool IntR_ModelBuildingNew(double *ei,double *Ii,int n,double rbPar[3]);
bool IntR_ModelBuildingNew(vector<double> &ei,vector<double> &Ii,double rbPar[3]);
double getpValue(double x,double mu,double a,double b,double ei);
double my_fit_linear(vector<double> &x,vector<double> &y,double c[2]);

double my_mean(double *x,size_t n);
double my_mean(vector<double> &x);
double my_STD(double *x,double m,size_t n);
double my_STD(vector<double> &x,double m);
double my_STD(vector<double> &x);
double my_min(double *x,size_t n);
double my_min(vector<double> &x);
double my_max(double *x,size_t n);
double my_max(vector<double> &x);
int my_count(double *x,double min_lim,double max_lim,size_t n);
int my_count(vector<double> &x,double min_lim,double max_lim);
bool my_RobustNormalFit(double *x,size_t n,double *mu,double *sigma,double *PW);
bool my_RobustNormalFit(vector<double> &x,double *mu,double *sigma,double *PW);
double my_sum(double *x,size_t n);
double my_sum(vector<double> &x);
double normpdf(double x,double mu,double sigma);
void normpdfV(double *x,size_t n,double *P,double mu,double sigma);
void normpdfV(vector<double> &x,vector<double> &P,double mu,double sigma);
