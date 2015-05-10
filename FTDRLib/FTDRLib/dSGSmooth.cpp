#include "StdAfx.h"
#include "dSGSmooth.h"
#include "OTrace.h"

dSGSmooth::dSGSmooth(void)
{
	X=NULL;
	dy=NULL;
	c=NULL;
	cov=NULL;
	work=NULL;

	datapoints=0;
	level=1;
	currentIDX=0;
	IsInitial=false;
	cenx=0;
}

void dSGSmooth::freeMEM()
{
	if(X!=NULL) gsl_matrix_free(X);
	if(dy!=NULL) gsl_vector_free(dy);
	if(c!=NULL) gsl_vector_free(c);
	if(cov!=NULL) gsl_matrix_free(cov);	
	if(work!=NULL) gsl_multifit_linear_free(work);	
}

dSGSmooth::~dSGSmooth(void)
{
	freeMEM();
}

//generally, n is a prime:2*j+1
bool dSGSmooth::Iinital(double *x,double *y,int n,int k)
{
	if(k<=0||n<=0) return false;
	if(k>=n) return false;
	level=k;
	datapoints=n;
	int p=level+1;	
	int i,j;
	freeMEM();
	cenx=x[datapoints/2+1];

	cov=gsl_matrix_alloc(p,p);
	X=gsl_matrix_alloc(n,p);
	c=gsl_vector_alloc(p);
	dy=gsl_vector_alloc(n);
	for(i=0;i<n;i++)
	{	
		double tx=1;
		double ttx=x[i]-cenx;
		for(j=0;j<=level;j++)
		{
			gsl_matrix_set(X,i,j,tx);
			tx*=ttx;
		}
		gsl_vector_set(dy,i,y[i]);
	}

	double chisq;
	work=gsl_multifit_linear_alloc(n,p);
	gsl_multifit_linear(X,dy,c,cov,&chisq,work);	
	currentIDX=0;
	IsInitial=true;
	return true;
}

bool dSGSmooth::Iinital(int n,int k)
{
	if(k<=0||n<=0) return false;
	if(k>=n) return false;
	level=k;
	datapoints=n;
	int p=level+1;		
	freeMEM();

	cov=gsl_matrix_alloc(p,p);
	X=gsl_matrix_alloc(n,p);
	gsl_matrix_set_zero(X);
	c=gsl_vector_alloc(p);
	dy=gsl_vector_alloc(n);
	work=gsl_multifit_linear_alloc(n,p);
	currentIDX=0;
	IsInitial=false;
	cenx=0;
	return true;
}

void dSGSmooth::ReInitial()
{
	currentIDX=0;
	IsInitial=false;
	cenx=0;
}

void dSGSmooth::Add(double x,double y)
{
	int i,j;
	double tx,ttx;
	j=(currentIDX+datapoints/2+2)%datapoints;
	double cennew=gsl_matrix_get(X,j,1)+cenx;
	for(i=0;i<datapoints;i++)
	{
		tx=gsl_matrix_get(X,i,1);
		tx+=cenx;
		tx-=cennew;
		ttx=tx;
		for(j=1;j<=level;j++) 
		{
			gsl_matrix_set(X,i,j,ttx);
			ttx*=tx;
		}
	}
	cenx=cennew;

	tx=1;
	ttx=x-cenx;
	for(int i=0;i<=level;i++)
	{
		gsl_matrix_set(X,currentIDX,i,tx);
		tx*=ttx;
	}
	gsl_vector_set(dy,currentIDX,y);
	currentIDX++;	
	if((!IsInitial)&&(currentIDX==datapoints))IsInitial=true;	
	double chisq;
	if(IsInitial) 
	{
		gsl_multifit_linear(X,dy,c,cov,&chisq,work);
		//printf("Update coff: ");
		/*for(int i=0;i<=level;i++)
		{
			printf("%lf\t",gsl_vector_get(c,i));
		}
		printf("\n");*/			
	}
	currentIDX=currentIDX%datapoints;	
}

double dSGSmooth::GetValue()
{
	if(!IsInitial) return -1;
	int idx=(currentIDX+datapoints/2+1)%datapoints;
	double br=gsl_vector_get(c,0);
	//for(int i=0;i<=level;i++)
	//{
	//	br+=gsl_vector_get(c,i)*gsl_matrix_get(X,idx,i);	
	//}
	//printf("compare: y=%lf\tsy=%lf\n",br,gsl_vector_get(dy,idx));		
	return br;
}

double dSGSmooth::GetCen()
{
	return cenx;
}


DiffCut::DiffCut()
{
	diffdata=NULL;
	xdata=NULL;
	currentIDX=0;	
	storeNum=0;
	dgt.Iinital(SW,LW);	
	IsInitial=false;
}

DiffCut::~DiffCut()
{
	if(diffdata!=NULL) delete []diffdata;
	if(xdata!=NULL) delete []xdata;
}

bool DiffCut::Initial(int n)
{
	if(n<=0) return false;
	if(diffdata!=NULL) delete []diffdata;
	if(xdata!=NULL) delete []xdata;
	diffdata=new double[n];
	xdata=new double[n];
	for(int i=0;i<n;i++)
	{
		diffdata[i]=0;
		xdata[i]=0;
	}
	storeNum=n;
	currentIDX=0;
	IsInitial=false;
	dgt.Iinital(SW,LW);	
	return true;
}

bool DiffCut::Initial(int n,int sw,int lw)
{
	if(n<=0||sw<lw||lw<1) return false;
	if(diffdata!=NULL) delete []diffdata;
	if(xdata!=NULL) delete []xdata;
	diffdata=new double[n];
	xdata=new double[n];
	for(int i=0;i<n;i++)
	{
		diffdata[i]=0;
		xdata[i]=0;
	}
	storeNum=n;
	currentIDX=0;
	IsInitial=false;
	dgt.Iinital(sw,lw);	
	return true;

}

void DiffCut::ReInitial()
{
	for(int i=0;i<storeNum;i++)
	{
		diffdata[i]=0;
		xdata[i]=0;
	}
	currentIDX=0;
	IsInitial=false;
	dgt.ReInitial();	
}

void DiffCut::Add(double x,double y)
{
	dgt.Add(x,y);
	double val=dgt.GetValue();
	if(val>0)
	{
		diffdata[currentIDX]=val;
		xdata[currentIDX]=dgt.GetCen();
		currentIDX++;
		if(!IsInitial)
		{
			if(currentIDX==storeNum) IsInitial=true;
		}
		currentIDX%=storeNum;
	}
}

bool DiffCut::IsMinmal()
{
	if(!IsInitial) return false;
	double *diffbuf;
	diffbuf=new double[storeNum];
	int i=currentIDX;
	int j=currentIDX+1;
	j%=storeNum;
	for(int k=0;k<storeNum-1;k++)
	{
		diffbuf[k]=diffdata[j]-diffdata[i];
		j++;
		i++;
		i=i%storeNum;
		j=j%storeNum;
	}
	
	i=storeNum/2;
	bool Isminp=true;
	for(j=0;j<i;j++)
	{
		if(diffbuf[j]>1e-6) 
		{
			Isminp=false;
			break;
		}
	}
	for(j=i;j<storeNum-1;j++) 
	{
		if(diffbuf[j]<-1e-6) 
		{
			Isminp=false;
			break;
		}
	}
	delete []diffbuf;
	return Isminp;
}

double DiffCut::GetMinmal()
{
	int i=(storeNum/2+1+currentIDX)%storeNum;
	return xdata[i];
}