#include "StdAfx.h"
#include "ISGSmooth.h"

ISGSmooth::ISGSmooth(void)
{
	xdata=NULL;
	ydata=NULL;
	currentIDX=0;
	c=NULL;
	val=NULL;
	vald=NULL;
	datapoints=0;
	level=0;
}

void ISGSmooth::freeMEM()
{
	if(xdata!=NULL) delete []xdata;
	if(ydata!=NULL) delete []ydata;
	if(c!=NULL) delete []c;
	if(val!=NULL) delete []val;
	if(vald!=NULL) delete []vald;
}

ISGSmooth::~ISGSmooth(void)
{
	freeMEM();
}

bool ISGSmooth::Initial(double *x,double *y,int n,int k)
{
	if(k<=0||n<=0) return false;
	if(k>=n) return false;
	level=k;
	datapoints=n;
	int p=level+1;	
	int i,j;
	freeMEM();

	val=gsl_matrix_alloc(p,p);
	c=gsl_vector_alloc(p);

	gsl_matrix *X;
	gsl_vector *dy;
	X=gsl_matrix_alloc(n,p);	
	dy=gsl_vector_alloc(n);	

	for(i=0;i<n;i++)
	{
		xdata[i]=x[i];
		ydata[i]=y[i];

		gsl_matrix_set(X,i,0,1);
		double tx=x[i];
		for(j=1;j<=level;j++)
		{
			gsl_matrix_set(X,i,j,tx);
			tx*=x[i];
		}
		gsl_vector_set(dy,i,y[i]);
	}
	double chisq;
	gsl_multifit_linear_workspace * work;
	work=gsl_multifit_linear_alloc(n,p);
	gsl_multifit_linear(X,dy,c,val,&chisq,work);
	gsl_multifit_linear_free(work);	
	gsl_matrix_free(X);
	gsl_vector_free(dy);
	gsl_matrix_scale(val,1/chisq);	
	gsl_matrix_memcpy(vald,val);
	currentIDX=0;
	IsInitial=true;
	return true;
}

void ISGSmooth::m_vector_mul(gsl_vector *x,gsl_vector *y,gsl_matrix *T)
{
	size_t sg=x->size;
	for(size_t i=0;i<sg;i++)
	{
		for(size_t k=0;k<sg;k++)
		{
			double mul=gsl_vector_get(x,i)*gsl_vector_get(y,k);
			gsl_matrix_set(T,i,k,mul);
		}
	}
}

bool ISGSmooth::add(double x,double y)
{
	if(!IsInitial) return false;
	//add new x
	gsl_vector *fi;
	gsl_vector *gn;
	gsl_vector *bt;
	gsl_matrix *dp;
	gsl_matrix *dtpn;
	int p=level+1;
	fi=gsl_vector_alloc(p);
	gn=gsl_vector_alloc(p);
	bt=gsl_vector_alloc(p);
	dp=gsl_matrix_alloc(p,p);
	dtpn=gsl_matrix_alloc(p,p);
	double tx=1;
	for(int i=0;i<p;i++)
	{
		gsl_vector_set(fi,i,tx);
		tx*=x;
	}
	gsl_blas_dgemv(CblasNoTrans,1.0,val,fi,0,gn);
	double beta;
	gsl_blas_ddot(fi,gn,&beta);
	beta=1/(1+beta);
	gsl_vector_scale(gn,beta);
	//update parameters
	gsl_blas_ddot(fi,c,&beta);
	beta=y-beta;
	gsl_vector_memcpy(bt,gn);
	gsl_vector_scale(bt,beta);
	gsl_vector_add(c,bt);
	///update val
	m_vector_mul(gn,fi,dp);
	gsl_blas_dsymm(CblasRight,CblasUpper,-1.0,val,dp,0,dtpn);
	gsl_matrix_add(val,dtpn);

	//remomove old x
	gsl_blas_dgemv(CblasNoTrans,1.0,vald,fi,0,gn);
	gsl_blas_ddot(fi,gn,&beta);
	beta=1/(1-beta);
	gsl_vector_scale(gn,beta);
	//
	gsl_blas_ddot(fi,c,&beta);
	beta=beta-ydata[currentIDX];
	gsl_vector_memcpy(bt,gn);
	gsl_vector_scale(bt,beta);
	gsl_vector_add(c,bt);
	//
	m_vector_mul(gn,fi,dp);
	gsl_blas_dsymm(CblasRight,CblasUpper,1.0,vald,dp,0,dtpn);
	gsl_matrix_add(vald,dtpn);
	xdata[currentIDX]=x;
	ydata[currentIDX]=y;
	currentIDX++;
	currentIDX=currentIDX%datapoints;
	return true;
}

double ISGSmooth::GetValue()
{
	int idx=currentIDX+datapoints/2+1;
	idx%=datapoints;
	double br=0;
	double tx=1;
	for(int i=0;i<=level;i++)
	{
		br+=gsl_vector_get(c,i)*tx;	
		tx*=xdata[idx];
	}
	return br;

}