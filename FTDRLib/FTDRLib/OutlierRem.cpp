#include "stdafx.h"
#include "OutlierRem.h"
#include "math.h"
//#define GSL_DLL
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

OutlierRem::OutlierRem(void)
{
	_K=2.698;
}

OutlierRem::~OutlierRem(void)
{
}

double OutlierRem::_mean()
{
	size_t size=_x.size();
	double m=0;
	size_t CT=0;
	for(size_t i=0;i<size;i++)
	{
		m+=_w[i]*_x[i];
		CT+=_w[i];
	}
	if(CT>0) m/=CT;
	return m;
}

double OutlierRem::_std(double m)
{
	size_t size=_x.size();	
	size_t CT=0;
	double std=0;
	for(size_t i=0;i<size;i++)
	{		
		std+=_w[i]*(_x[i]-m)*(_x[i]-m);
		CT+=_w[i];
	}
	if(CT>0) std/=CT;
	return sqrt(std);
}

double OutlierRem::Mean(vector<double> &x)
{
	size_t size=x.size();
	double m=0;
	for(size_t i=0;i<size;i++)m+=x[i];	
	if(size>0) m/=size;
	return m;
}

double OutlierRem::Std(double m,vector<double> &x)
{
	size_t size=x.size();	
	double std=0;
	for(size_t i=0;i<size;i++)std+=(x[i]-m)*(x[i]-m);
	if(size>0) std/=size;
	return sqrt(std);
}

double OutlierRem::MeanW(vector<double> &x,vector<double> &w)
{
	size_t size=x.size();
	if(size<=0) return 0;
	double m=0;
	double sumw=0;
	double mb=0;
	for(size_t i=0;i<size;i++)
	{
		m+=x[i]*w[i];
		mb+=x[i];
		sumw+=w[i];
	}
	if(sumw>1e-6) m/=sumw;
	else m=mb/size;
	return m;
}

double OutlierRem::StdW(double m,vector<double> &x,vector<double> &w)
{
	size_t size=x.size();
	if(size<=0) return 0;
	double std=0;
	double sumw=0;
	double stdb=0;
	for(size_t i=0;i<size;i++)
	{
		std+=(x[i]-m)*(x[i]-m)*w[i];
		stdb+=(x[i]-m)*(x[i]-m);
		sumw+=w[i];
	}
	if(sumw>1e-6) std/=sumw;
	else std=stdb/size;
	return sqrt(std);
}

bool OutlierRem::_updatew()
{
	double m=_mean();
	double std=_std(m);
	size_t n=_x.size();
	double cutoff=_K*std;
	for(size_t i=0;i<n;i++)
	{
		if(fabs(_x[i]-m)>cutoff) _w[i]=0;
		else _w[i]=1;
	}
	m_new=_mean();
	std_new=_std(m_new);
	double erro=abs(m_new-m)+abs(std_new-std);
	if(erro<1e-4) return true;
	return false;
}

void OutlierRem::initialw(vector<double> *data)
{
	_w.clear();
	_x.clear();
	size_t size=data->size();
	for(size_t i=0;i<size;i++)
	{
		_w.push_back(1);
		_x.push_back(data->at(i));
	}
}

bool OutlierRem::oRemove(vector<double> &x)
{
	initialw(&x);
	size_t step=0;
	while(step<100)
	{
		if(_updatew()) break;
		step++;
	}
	if(step<100)//sucessful
	{
		size_t n=_x.size();
		x.clear();
		for(size_t i=0;i<n;i++)
			if(_w[i]) x.push_back(_x[i]);
		return true;
	}	
	return false;
}

double OutlierRem::GetMean()
{
	return m_new;
}

double OutlierRem::GetStd()
{
	return std_new;
}

//Grubbs' test for outliers:http://en.wikipedia.org/wiki/Grubbs'_test_for_outliers
bool OutlierRem::oRemoveGTest(vector<double> &x)
{	
	size_t N=x.size();
	if(N<=4) return false;
	double *y,*p;
	y=new double[N];
	int MaximalTry=N/4;
	int i,n;
	for(i=0;i<N;i++) y[i]=x[i];
	gsl_sort(y,1,N);
	size_t step=0;
	int begin=0;
	int end=N-1;	
	bool Is_C=true;
	while(step<MaximalTry&&Is_C)
	{
		p=y+begin;
		n=end-begin+1;	
		m_new=gsl_stats_mean(p,1,n);
		std_new=gsl_stats_sd_m(p,1,n,m_new);
		if(std_new<=1e-6) break;
		double g_min=(m_new-y[begin])/std_new;
		double g_max=(y[end]-m_new)/std_new;
		double ta=0.025/N;
		ta=gsl_cdf_tdist_Qinv(ta,n-2);	
		ta=ta*ta;
		ta=(n-1)*sqrt(ta/(n*(n-2+ta)));
		Is_C=false;
		if(g_min>ta) 
		{
			begin++;
			Is_C=true;
		}
		if(g_max>ta) 
		{
			end--;
			Is_C=true;
		}
		step++;		
	}	
	x.clear();
	for(i=begin;i<=end;i++) x.push_back(y[i]);
	delete []y;
	return true;
}


bool OutlierRem::oRemoveGTestW(vector<double> &x,vector<double> &w)
{
	size_t N=x.size();
	if(N<=4) return false;
	double *y,*p;
	y=new double[N];
	int MaximalTry=N/4;
	int i,n;	
	for(i=0;i<N;i++) y[i]=x[i];
	size_t *ID;
	ID=new size_t[N];
	gsl_sort_index(ID,y,1,N);

	size_t step=0;
	int begin=0;
	int end=N-1;	
	bool Is_C=true;
	while(step<MaximalTry&&Is_C)
	{
		p=y+begin;
		n=end-begin+1;	
		m_new=gsl_stats_mean(p,1,n);
		std_new=gsl_stats_sd_m(p,1,n,m_new);
		if(std_new<=1e-6) break;
		double g_min=(m_new-y[begin])/std_new;
		double g_max=(y[end]-m_new)/std_new;
		double ta=0.025/N;
		ta=gsl_cdf_tdist_Qinv(ta,n-2);	
		ta=ta*ta;
		ta=(n-1)*sqrt(ta/(n*(n-2+ta)));
		Is_C=false;
		if(g_min>ta) 
		{
			begin++;
			Is_C=true;
		}
		if(g_max>ta) 
		{
			end--;
			Is_C=true;
		}
		step++;		
	}
	
	x.clear();
	m_new=0;
	double sumw=0;
	for(i=begin;i<=end;i++) 
	{
		x.push_back(y[i]);
		m_new+=y[i]*w[ID[i]];
		sumw+=w[ID[i]];
	}
	if(sumw>1e-4) m_new/=sumw;
	std_new=0;
	for(i=begin;i<=end;i++) 
	{		
		std_new+=(y[i]-m_new)*(y[i]-m_new)*w[ID[i]];		
	}
	if(sumw>1e-4) std_new=sqrt(std_new/sumw);
	//m_new=gsl_stats_mean(p,1,n);
	//std_new=gsl_stats_sd_m(p,1,n,m_new);	
	delete []y;
	delete []ID;
	return true;
}

//Grubbs' test for outliers:http://en.wikipedia.org/wiki/Grubbs'_test_for_outliers
bool OutlierRem::oRemoveGTest(double *y, int N)
{		
	if(N<=4) return false;
	double *p;
	int MaximalTry=N/4;
	int i,n;	
	gsl_sort(y,1,N);
	size_t step=0;
	int begin=0;
	int end=N-1;	
	bool Is_C=true;
	while(step<MaximalTry&&Is_C)
	{
		p=y+begin;
		n=end-begin+1;	
		m_new=gsl_stats_mean(p,1,n);
		std_new=gsl_stats_sd_m(p,1,n,m_new);
		if(std_new<=1e-6) break;
		double g_min=(m_new-y[begin])/std_new;
		double g_max=(y[end]-m_new)/std_new;
		double ta=0.025/N;
		ta=gsl_cdf_tdist_Qinv(ta,n-2);	
		ta=ta*ta;
		ta=(n-1)*sqrt(ta/(n*(n-2+ta)));
		Is_C=false;
		if(g_min>ta) 
		{
			begin++;
			Is_C=true;
		}
		if(g_max>ta) 
		{
			end--;
			Is_C=true;
		}
		step++;		
	}
	if(begin==0) return true;
	for(i=0;i<n;i++)  y[i]=y[i+begin];//move back	
	return true;
}

//Chauvenet's criterion:  http://en.wikipedia.org/wiki/Chauvenet%27s_criterion
bool OutlierRem::oRemoveChauCt(vector<double> &x)
{	
	size_t N=x.size();
	if(N<=4) return false;
	double *y,*p;
	y=new double[N];
	int MaximalTry=N/4;
	int i,n;
	for(i=0;i<N;i++) y[i]=x[i];
	gsl_sort(y,1,N);
	size_t step=0;
	int begin=0;
	int end=N-1;	
	while(step<MaximalTry)
	{
		p=y+begin;
		n=end-begin+1;	
		m_new=gsl_stats_mean(p,1,n);
		std_new=gsl_stats_sd_m(p,1,n,m_new);
		if(std_new<=1e-6) break;
		double CutOff=0.5/n;
		for(;begin<end;begin++)
		{
			double z=fabs(y[begin]-m_new);
			z=gsl_cdf_gaussian_Q(z,std_new);
			if(z>=CutOff) break;
		}		
		for(;end>begin;end--)
		{
			double z=fabs(y[begin]-m_new);
			z=gsl_cdf_gaussian_Q(z,std_new);
			if(z>=CutOff) break;
		}
		if(end-begin+1==n) break;
		step++;		
	}	
	x.clear();
	for(i=begin;i<=end;i++) x.push_back(y[i]);
	delete []y;
	return true;
}

bool OutlierRem::oRemoveChauCt(double *y, int N)
{	
	if(N<=4) return false;
	double *p;
	int MaximalTry=N/4;
	int i,n;
	gsl_sort(y,1,N);
	size_t step=0;
	int begin=0;
	int end=N-1;	
	bool Is_C=true;
	while((step<MaximalTry)&&Is_C)
	{
		p=y+begin;
		n=end-begin+1;	
		m_new=gsl_stats_mean(p,1,n);
		std_new=gsl_stats_sd_m(p,1,n,m_new);
		if(std_new<=1e-6) break;
		double CutOff=0.5/n;
		for(;begin<end;begin++)
		{
			double z=fabs(y[begin]-m_new);
			z=gsl_cdf_gaussian_Q(z,std_new);
			if(z>=CutOff) break;
		}		
		for(;end>begin;end--)
		{
			double z=fabs(y[begin]-m_new);
			z=gsl_cdf_gaussian_Q(z,std_new);
			if(z>=CutOff) break;
		}
		if(end-begin+1==n) break;
		step++;		
	}	
	if(begin==0) return true;
	for(i=0;i<n;i++)  y[i]=y[i+begin];//move back
	return true;
}

int OutlierRem::count(double *x,int n,double min_v, double max_v)
{
	int c=0;
	for(int i=0;i<n;i++)
	{
		if(x[i]>min_v&&x[i]<max_v) c++;
	}
	return c;
}

double OutlierRem::gsl_sum(double *x,int n)
{
	double sumv=0;
	for(int i=0;i<n;i++)sumv+=x[i];
	return sumv;
}

bool OutlierRem::ORemMixModel(double *x,int n)
{	
	double mu=gsl_stats_mean(x,1,n);
	double sigma=gsl_stats_sd_m(x,1,n,mu);
	double a=gsl_stats_min(x,1,n);
	double b=gsl_stats_max(x,1,n);
	double alim=mu-2*(sigma);
	double blim=mu+2*(sigma);
	double L=alim-a+b-blim;
	double K=b-a;
	if(L<=0) return false; 
	size_t deltn=n-count(x,n,alim,blim);	
	double PW2=K*deltn/(L*n);
	double PW1=1-PW2;
	double muold=mu;
	double sigmaold=sigma;
	double e=1;
	int iter=0;
	double *PWP1,*PWAll;
	PWP1=new double[n];
	PWAll=new double[n];
	size_t i;
	while(e>0.00001)  
	{
		for(i=0;i<n;i++)
		{
			PWP1[i]=PW1*gsl_ran_gaussian_pdf(x[i]-muold,sigmaold);
			PWAll[i]=PWP1[i]+PW2/K;
			PWP1[i]=PWP1[i]/PWAll[i];
		}		
		PW1=gsl_sum(PWP1,n);
		mu=0;
		for(i=0;i<n;i++)mu+=PWP1[i]*x[i];		
		mu/=PW1;
		sigma=0;
		for(i=0;i<n;i++)sigma+=PWP1[i]*(x[i]-muold)*(x[i]-muold);
		sigma=sqrt(sigma/PW1);
		PW1=PW1/n;
		PW2=1-PW1;
		e=fabs(sigma-sigmaold)+fabs(mu-muold);
		muold=mu;
		sigmaold=sigma;
		iter=iter+1;
		if(iter>200) break;
	}
	delete []PWP1;
	delete []PWAll;
	m_new=mu;
	std_new=sigma;
	return true;
}

bool OutlierRem::ORemMixModel(vector<double> &x)
{
	size_t n=x.size();
	if(n<=0) return false;
	else if(n<=3)//special case, directly std and mean
	{
		m_new=Mean(x);
		std_new=Std(m_new,x);
		return true;
	}
	double *y;
	y=new double[n];
	for(int i=0;i<n;i++) y[i]=x[i];
	bool bRt=ORemMixModel(y,n);
	delete []y;
	return bRt;
}


