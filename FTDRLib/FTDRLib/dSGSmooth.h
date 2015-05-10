#pragma once

#include "GaussFit.h"
#define SW 21
#define LW 2

class dSGSmooth
{
private:
	gsl_multifit_linear_workspace * work;
	gsl_matrix *X;
	gsl_vector *dy;
	gsl_vector *c;
	gsl_matrix *cov;	
	int datapoints;
	int level;
	int currentIDX;
	void freeMEM();
	bool IsInitial;
	double cenx;
public:
	dSGSmooth(void);
	~dSGSmooth(void);
	bool Iinital(double *x,double *y,int n,int k);
	bool Iinital(int n,int k);
	void ReInitial();
	void Add(double x,double y);
	double GetValue();
	double GetCen();
};


class DiffCut
{
private:
	double *diffdata;
	double *xdata;
	int currentIDX;	
	int storeNum;
	dSGSmooth dgt;
	bool IsInitial;
public:
	DiffCut();
	~DiffCut();
	bool Initial(int n);
	bool Initial(int n,int sw,int lw);
	void ReInitial();
	void Add(double x,double y);
	bool IsMinmal();
	double GetMinmal();
};