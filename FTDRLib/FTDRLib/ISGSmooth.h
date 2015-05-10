#pragma once
#include "GaussFit.h"

class ISGSmooth
{
private:
	double *xdata;
	double *ydata;
	int currentIDX;
	gsl_vector *c;
	gsl_matrix *val;
	gsl_matrix *vald;
	int datapoints;
	int level;
	void freeMEM();
	bool IsInitial;
	void m_vector_mul(gsl_vector *x,gsl_vector *y,gsl_matrix *T);
public:
	ISGSmooth(void);
	~ISGSmooth(void);
	bool Initial(double *x,double *y,int n,int k);
	bool add(double x,double y);
	double GetValue();
};
