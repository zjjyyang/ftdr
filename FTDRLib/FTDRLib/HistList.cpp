#include "stdafx.h"
#include "HistList.h"
#include "gsl\gsl_statistics_double.h"
#include "math.h"

myHistogram::myHistogram()
{
	range=NULL;
	bin=NULL;
	n=0;
	mean=0;
	std=0;
}

myHistogram::myHistogram(const myHistogram &r)
{
	n=r.n;
	if(n>0)
	{
		range=new double[n+1];
		bin=new double[n];
		for(int i=0;i<n;i++) 
		{
			range[i]=r.range[i];
			bin[i]=r.bin[i];
		}
		range[n]=r.range[n];
	}
	mean=r.mean;
	std=r.std;
	rawname=r.rawname;
}

myHistogram::myHistogram(gsl_histogram *r)
{
	if(n>0) 
	{
		if(range!=NULL) delete []range;
		if(bin!=NULL) delete []bin;
		range=NULL;
		bin=NULL;
	}
	n=r->n;
	if(n>0)
	{
		range=new double[n+1];
		bin=new double[n];
		for(int i=0;i<n;i++) 
		{
			range[i]=r->range[i];
			bin[i]=r->bin[i];
		}
		range[n]=r->range[n];
	}		
}

myHistogram::~myHistogram()
{
	if(n>0) 
	{
		if(range!=NULL) delete []range;
		if(bin!=NULL) delete []bin;
		range=NULL;
		bin=NULL;
		n=0;
	}
}

myHistogram &myHistogram::operator=(const myHistogram &r)
{
	if(n>0) 
	{
		if(range!=NULL) delete []range;
		if(bin!=NULL) delete []bin;
		range=NULL;
		bin=NULL;		
	}
	n=r.n;
	if(n>0)
	{
		range=new double[n+1];
		bin=new double[n];
		for(int i=0;i<n;i++) 
		{
			range[i]=r.range[i];
			bin[i]=r.bin[i];
		}
		range[n]=r.range[n];
	}
	mean=r.mean;
	std=r.std;
	rawname=r.rawname;
	return *this;
}

double myHistogram::GetMaxY()
{
	if(n<=0) return 0;
	double maxy=0;
	for(int i=0;i<n;i++)
	{
		if(maxy<bin[i]) maxy=bin[i];
	}
	if(maxy<C0/std) maxy=C0/std;
	return maxy;
}

double myHistogram::GetMinX()
{
	if(n>0)	return range[0];
	return -1;
}

double myHistogram::GetMaxX()
{
	if(n>0)	return range[n];
	return -1;
}

double myHistogram::GetProbGS(double x)
{
	if(std>1e-6)
	{
		double t=(x-mean)/std;
		t=t*t/2;
		return C0*exp(-t)/std;
	}
	return 0;
}

bool myHistogram::alloc(int num)
{
	if(num<=0) return false;
	if(n>0) 
	{
		if(range!=NULL) delete []range;
		if(bin!=NULL) delete []bin;
		range=NULL;
		bin=NULL;
		n=0;
	}
	n=num;
	range=new double[n+1];
	bin=new double[n];
	return true;
}

bool myHistogram::Convert(double *res,int n)
{
	gsl_histogram *r;
	if(n<=2) return false;
	double std=gsl_stats_sd(res,1,n);
	double bin=3.49*std/pow(n*1.0,1.0/3);//Scott's ruler
	if(bin<=0) return false;	
	double a=gsl_stats_min(res,1,n);
	double b=gsl_stats_max(res,1,n);
	int num=(int)((b-a)/bin);
	r=gsl_histogram_alloc(num);
	gsl_histogram_set_ranges_uniform(r,a,b);
	for(int i=0;i<n;i++)
	{
		gsl_histogram_increment(r,res[i]);
	}	
	Convert(r,n);	
	gsl_histogram_free(r);	
	return true;
}

bool myHistogram::Convert(vector<double> &x)
{
	gsl_histogram *r;
	size_t n=x.size();
	if(n<=2) return false;
	double *res;
	res=new double[n];
	size_t i;
	for(i=0;i<n;i++) res[i]=x[i];
	double std=gsl_stats_sd(res,1,n);
	double bin=3.49*std/pow(n*1.0,1.0/3);//Scott's ruler
	if(bin<=0) 
	{
		delete []res;
		return false;
	}
	double a=gsl_stats_min(res,1,n);
	double b=gsl_stats_max(res,1,n);
	int num=(int)((b-a)/bin);
	r=gsl_histogram_alloc(num);
	gsl_histogram_set_ranges_uniform(r,a,b);
	for(i=0;i<n;i++)
	{
		gsl_histogram_increment(r,res[i]);
	}	
	Convert(r,n);	
	gsl_histogram_free(r);
	delete []res;	
	return true;
}

void myHistogram::Convert(gsl_histogram *r,int num)
{
	if(num<=0) return;
	if(n>0) 
	{
		if(range!=NULL) delete []range;
		if(bin!=NULL) delete []bin;
		range=NULL;
		bin=NULL;
		n=0;
	}
	n=r->n;
	if(n>0)
	{
		double w=r->range[1]-r->range[0];
		range=new double[n+1];
		bin=new double[n];
		for(int i=0;i<n;i++) 
		{
			range[i]=r->range[i];
			bin[i]=r->bin[i]/(num*w);
		}
		range[n]=r->range[n];//n+1
	}
}

HistList::HistList(void)
{	
	Acount=PRE_MC;
	HPL=new myHistogram[PRE_MC];
	Mcount=0;
}

HistList::HistList(int n)
{
	Mcount=0;
	if(n>0)
	{
		HPL= new myHistogram[n];
		Acount=n;
	}
	else
	{	
		Acount=PRE_MC;
		HPL=new myHistogram[PRE_MC];
	}
}

HistList::~HistList(void)
{
	if(Acount>0&&HPL!=NULL) delete []HPL;
}

bool HistList::AddHist(gsl_histogram *r,double m,double sd,CString raw)
{
	if(r==NULL) return false;	
	int n=r->n;
	if(n<=0) return false;
	myHistogram mt(r);	
	mt.mean=m;
	mt.std=sd;
	mt.rawname=raw;
	push_back(mt);
	return true;
}

int HistList::size()
{
	return Mcount;
}

bool HistList::at(int pos,myHistogram &mt)
{
	if(pos<0||pos>Mcount) return false;
	mt=HPL[pos];
	return true;
}

HistList::HistList(HistList &r)
{
	if(Acount>0&&HPL!=NULL)
	{
		delete []HPL;
		HPL=NULL;
		Acount=0;
		Mcount=0;
	}	
	int count=r.size();
	myHistogram mt;
	for(int i=0;i<count;i++)
	{
		if(r.at(i,mt))push_back(mt);
	}
}

void HistList::push_back(const myHistogram mt)
{
	int i;
	if(Acount<=Mcount)
	{
		myHistogram *tpl;
		Acount+=ADD_CT;
		tpl=new myHistogram[Acount];
		for(i=0;i<Mcount;i++)
		{
			tpl[i]=HPL[i];
		}
		tpl[Mcount]=mt;
		Mcount++;
		if(Mcount>1) delete []HPL;
		HPL=tpl;
	}
	else
	{
		HPL[Mcount]=mt;
		Mcount++;
	}
}


HistList &HistList::operator=(HistList &r)
{
	if(Acount>0&&HPL!=NULL)
	{
		delete []HPL;
		HPL=NULL;
		Acount=0;
		Mcount=0;
	}	
	int count=r.size();
	myHistogram mt;
	for(int i=0;i<count;i++)
	{
		if(r.at(i,mt))push_back(mt);
	}
	return *this;
}

bool HistList::AddHist(double *res,int n,double m,double sd,CString raw)
{
	gsl_histogram *r;
	if(n<=0) return false;
	double std=gsl_stats_sd(res,1,n);
	double bin=3.49*std/pow(n*1.0,1.0/3);//Scott's ruler
	if(bin<=0) return false;	
	double a=gsl_stats_min(res,1,n);
	double b=gsl_stats_max(res,1,n);
	int num=(int)((b-a)/bin);
	r=gsl_histogram_alloc(num);
	gsl_histogram_set_ranges_uniform(r,a,b);
	for(int i=0;i<n;i++)
	{
		gsl_histogram_increment(r,res[i]);
	}
	bool bResult=AddHist(r,m,sd,raw);
	gsl_histogram_free(r);
	return bResult;
}

void HistList::reshape()
{
	if(Acount==Mcount) return;
	if(Mcount<=0) return;	
	myHistogram *tpl;
	tpl=new myHistogram[Mcount];
	for(int i=0;i<Mcount;i++)
	{
		tpl[i]=HPL[i];
	}
	Acount=Mcount;
	delete []HPL;
	HPL=tpl;
}

void HistList::clear()
{
	Mcount=0;
}