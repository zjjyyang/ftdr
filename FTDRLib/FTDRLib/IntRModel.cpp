#include "StdAfx.h"
#include "IntRModel.h"
#include "OTrace.h"
#include "OutlierRem.h"
#define PIS 0.398942280401433


FittinfPar::FittinfPar()
{
	a=0;
	b=0;
	theta=0;
}

FittinfPar::~FittinfPar()
{

}

FittinfPar::FittinfPar(const FittinfPar & r)
{
   size_t CT=r.yi.size();
   for(size_t i=0;i<CT;i++)
   {
	   yi.push_back(r.yi.at(i));
	   epsi.push_back(r.epsi.at(i));
	   weight.push_back(r.weight.at(i));
   }
   a=r.a;
   b=r.b;
   theta=r.theta;
}

FittinfPar &FittinfPar::operator =(const FittinfPar &r)
{
	epsi.clear();
	yi.clear();
	weight.clear();
   size_t CT=r.yi.size();
   for(size_t i=0;i<CT;i++)
   {
	   yi.push_back(r.yi.at(i));
	   epsi.push_back(r.epsi.at(i));
	   weight.push_back(r.weight.at(i));
   }
   a=r.a;
   b=r.b;
   theta=r.theta;
   return *this;
}

void FittinfPar::Initial()
{
	size_t sg=yi.size();
	for(size_t i=0;i<sg;i++) weight.push_back(1);
}

double FittinfPar::Error_mean()
{
	size_t size=yi.size();
	double m=0;
	int wt=0;
	for(size_t i=0;i<size;i++)
	{
		if(weight[i])
		{
			m+=epsi.at(i);
			wt++;
		}
	}
	if(wt>0) m/=wt;
	return m;
}

void FittinfPar::UpdateW()
{	
	double mu=Error_mean();
	size_t sg=yi.size();
	for(size_t i=0;i<sg;i++)
	{
		double sigma=a/sqrt(yi[i])+b;	
		double LB=mu-W_TUNE*sigma;
		double HB=mu+W_TUNE*sigma;
		if(epsi[i]>HB||epsi[i]<LB) weight[i]=0;
		else weight[i]=1;
	}
}

double FittinfPar::ParError(double aold,double bold)
{
	double as=fabs(a)+fabs(aold);
	double bs=fabs(b)+fabs(bold);
	double err=0;
	if(as>1e-4)err+=2*fabs(a-aold)/as;
	if(bs>1e-4) err+=2*fabs(b-bold)/bs;
	return err;
}

double NormalPDF(double x,double mu,double sigma)
{
	if(sigma<=0) return 0;
	double t=(x-mu)/sigma;
	return SQRTRPI*exp(-t*t/2)/sigma;
}

double LogNormalPDF(double x,double mu,double sigma)
{
	if(sigma<=1e-66) return 0;
	double t=(x-mu)/sigma;
	return LOG2PIH-log(sigma)-t*t/2;
}

double Error_mean(FittinfPar *Par)
{
	size_t size=Par->yi.size();
	double m=0;
	int wt=0;
	for(size_t i=0;i<size;i++)
	{
		if(Par->weight[i])
		{
			m+=Par->epsi.at(i);
			wt++;
		}
	}
	if(wt>0) m/=wt;
	return m;
}

double Error_std(FittinfPar *Par)
{
	size_t size=Par->yi.size();
	double m=Error_mean(Par);
	double std=0;
	int wt=0;
	for(size_t i=0;i<size;i++)
	{
		if(Par->weight[i])
		{
			double t=Par->epsi.at(i)-m;
			std+=t*t;
			wt++;
		}
	}
	if(wt>0) std/=wt;
	return std;
}

/*注意：为了满足求函数极小值的需要，这里求出的是-L，也就是极大似然的负值*/
//
//double IntR_Res_f(const gsl_vector *v, void *params)
//{
//	FittinfPar *Par;
//	Par=(FittinfPar *)params;
//	size_t count=Par->yi.size();
//	double mu=Error_mean(Par);
//	double sigma;
//	double L=0;
//	double a=gsl_vector_get(v,0);
//	double b=gsl_vector_get(v,1);
//	double theta=gsl_vector_get(v,2);
//	for(size_t i=0;i<count;i++)
//	{
//		double yi=Par->yi.at(i);
//		double x=Par->epsi.at(i);
//		sigma=a*exp(theta*yi)+b;		
//		L+=LogNormalPDF(x,mu,sigma);		
//	}	
//	return -L;
//}
//
//void IntR_Res_df(const gsl_vector *v, void *params,gsl_vector *df)
//{		
//	FittinfPar *Par;
//	Par=(FittinfPar *)params;
//	size_t count=Par->yi.size();
//	double mu=Error_mean(Par);
//	double sigma;
//	double L=0;
//	double a=gsl_vector_get(v,0);
//	double b=gsl_vector_get(v,1);
//	double theta=gsl_vector_get(v,2);
//
//	double DFa=0;
//	double DFb=0;
//	double DFTheta=0;
//	for(size_t i=0;i<count;i++)
//	{
//		double xi=Par->epsi.at(i);
//		double yi=Par->yi.at(i);		
//		sigma=a*exp(theta*yi)+b;	
//		double p2=1/sigma;
//		double p1=(xi-mu)*(xi-mu)*p2*p2*p2;
//		double p3=exp(theta*yi);
//		double ti=p1-p2;
//		DFb+=ti;
//		ti*=p3;
//		DFa+=ti;
//		DFTheta+=ti*yi;				
//	}	
//	gsl_vector_set(df,0,-DFa);//original is negtaive,-DFa
//	gsl_vector_set(df,1,-DFb);//-DFb
//	gsl_vector_set(df,2,-DFTheta*a);	//-DFTheta
//}	
//
//void IntR_Res_fdf(const gsl_vector *v, void *params,double *f,gsl_vector *df)
//{	
//	*f=IntR_Res_f(v,params);
//	IntR_Res_df(v,params,df);
//}

//新的模型，使用了MS次数作为变化量
//
double IntR_Res_f(const gsl_vector *v, void *params)
{
	FittinfPar *Par;
	Par=(FittinfPar *)params;
	size_t count=Par->yi.size();
	double mu=Error_mean(Par);
	double sigma;
	double L=0;
	double a=gsl_vector_get(v,0);
	double b=gsl_vector_get(v,1);	
	for(size_t i=0;i<count;i++)
	{
		if(Par->weight[i])
		{
			double yi=Par->yi.at(i);
			double x=Par->epsi.at(i);
			sigma=a/sqrt(yi)+b;		
			L+=LogNormalPDF(x,mu,sigma);	
		}
	}	
	return -L;
}

void IntR_Res_df(const gsl_vector *v, void *params,gsl_vector *df)
{		
	FittinfPar *Par;
	Par=(FittinfPar *)params;
	size_t count=Par->yi.size();
	double mu=Error_mean(Par);
	double sigma;
	double L=0;
	double a=gsl_vector_get(v,0);
	double b=gsl_vector_get(v,1);

	double DFa=0;
	double DFb=0;
	for(size_t i=0;i<count;i++)
	{
		if(Par->weight[i])
		{
			double xi=Par->epsi.at(i);
			double yi=Par->yi.at(i);		
			sigma=a/sqrt(yi)+b;	
			if(sigma>1e-16)
			{
				double p2=1/sigma;
				double p1=(xi-mu)*(xi-mu)*p2*p2*p2;
				double p3=1/sqrt(yi);
				double ti=p1-p2;
				DFb+=ti;
				ti*=p3;
				DFa+=ti;	
			}
		}
	}	
	gsl_vector_set(df,0,-DFa);//original is negtaive,-DFa
	gsl_vector_set(df,1,-DFb);//-DFb	
}	

void IntR_Res_fdf(const gsl_vector *v, void *params,double *f,gsl_vector *df)
{	
	*f=IntR_Res_f(v,params);
	IntR_Res_df(v,params,df);
}
//please make sure the Insnt[0]>eps before fitting
//double Fitting(FittinfPar *Par)
//{
//	size_t iter = 0;
//	int status;
//	const gsl_multimin_fdfminimizer_type *T;
//	gsl_multimin_fdfminimizer *s;
//	/* Starting point */
//	double a0=Error_std(Par)/3;
//	double b0=a0;
//	double theta0=-2;//-0.5;
//	a0*=30;
//	gsl_vector *v;
//	v=gsl_vector_alloc(3);
//	gsl_vector_set(v,0,a0);
//	gsl_vector_set(v,1,b0);
//	gsl_vector_set(v,2,theta0);
//	Par->a=a0;
//	Par->b=b0;
//	Par->theta=theta0;
//
//	gsl_multimin_function_fdf my_func;
//	my_func.f = &IntR_Res_f;
//	my_func.df = &IntR_Res_df;
//	my_func.fdf = &IntR_Res_fdf;
//	my_func.n = 3;
//	my_func.params =(void *)Par;
//
//	T = gsl_multimin_fdfminimizer_vector_bfgs;
//	s = gsl_multimin_fdfminimizer_alloc(T,3);//3 是未知数个数
//	gsl_multimin_fdfminimizer_set(s,&my_func,v,0.01,1e-5);
//	do
//	{
//		iter++;
//		status = gsl_multimin_fdfminimizer_iterate(s);
//		if(status==GSL_ENOPROG)
//		{	
//			Par->a=gsl_vector_get(s->x,0);
//			Par->b=gsl_vector_get(s->x,1);
//			Par->theta=gsl_vector_get(s->x,2);				
//			break;
//		}
//		status = gsl_multimin_test_gradient(s->gradient,1e-6);
//		if(status == GSL_SUCCESS)
//		{
//			Par->a=gsl_vector_get(s->x,0);
//			Par->b=gsl_vector_get(s->x,1);
//			Par->theta=gsl_vector_get(s->x,2);					
//		}
//	}while(status==GSL_CONTINUE&&iter<100);
//	//printf("Iter=%d\n",iter);
//	double breturn=IntR_Res_f(s->x,(void *)Par);
//	gsl_multimin_fdfminimizer_free(s);	
//	gsl_vector_free(v);
//	//cal the goodness of fitting
//	return -breturn;	
//}

double Fitting(FittinfPar *Par)
{
	size_t iter = 0;
	int status;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	/* Starting point */
	//adding on 2011.8.27
	Par->Initial();
	//////////////////
	double a0=Error_std(Par)*0.6;
	double b0=a0*0.4;
	gsl_vector *v;
	v=gsl_vector_alloc(2);
	gsl_vector_set(v,0,a0);
	gsl_vector_set(v,1,b0);
	Par->a=a0;
	Par->b=b0;	

	gsl_multimin_function_fdf my_func;
	my_func.f = &IntR_Res_f;
	my_func.df = &IntR_Res_df;
	my_func.fdf = &IntR_Res_fdf;
	my_func.n = 2;
	my_func.params =(void *)Par;

	T = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc(T,2);//3 是未知数个数
	gsl_multimin_fdfminimizer_set(s,&my_func,v,0.01,1e-5);
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
		if(status==GSL_ENOPROG)
		{	
			Par->a=gsl_vector_get(s->x,0);
			Par->b=gsl_vector_get(s->x,1);						
			break;
		}
		status = gsl_multimin_test_gradient(s->gradient,1e-6);
		if(status == GSL_SUCCESS)
		{
			Par->a=gsl_vector_get(s->x,0);
			Par->b=gsl_vector_get(s->x,1);							
		}
	}while(status==GSL_CONTINUE&&iter<100);
	//printf("Iter=%d\n",iter);
	double breturn=IntR_Res_f(s->x,(void *)Par);
	gsl_multimin_fdfminimizer_free(s);	
	gsl_vector_free(v);
	//cal the goodness of fitting
	return -breturn;	
}

double Fitting_W(FittinfPar *Par)
{
	size_t iter = 0;
	int status;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	/* Starting point */
	Par->UpdateW();
	gsl_vector *v;
	v=gsl_vector_alloc(2);
	gsl_vector_set(v,0,Par->a);
	gsl_vector_set(v,1,Par->b);
	//Par->a=a0;
	//Par->b=b0;	

	gsl_multimin_function_fdf my_func;
	my_func.f = &IntR_Res_f;
	my_func.df = &IntR_Res_df;
	my_func.fdf = &IntR_Res_fdf;
	my_func.n = 2;
	my_func.params =(void *)Par;

	T = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc(T,2);//3 是未知数个数
	gsl_multimin_fdfminimizer_set(s,&my_func,v,0.01,1e-5);
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
		if(status==GSL_ENOPROG)
		{	
			Par->a=gsl_vector_get(s->x,0);
			Par->b=gsl_vector_get(s->x,1);						
			break;
		}
		status = gsl_multimin_test_gradient(s->gradient,1e-6);
		if(status == GSL_SUCCESS)
		{
			Par->a=gsl_vector_get(s->x,0);
			Par->b=gsl_vector_get(s->x,1);							
		}
	}while(status==GSL_CONTINUE&&iter<100);
	//printf("Iter=%d\n",iter);
	double breturn=IntR_Res_f(s->x,(void *)Par);
	gsl_multimin_fdfminimizer_free(s);	
	gsl_vector_free(v);
	//cal the goodness of fitting
	return -breturn;	
}

bool IntR_ModelBuilding(double *ei,double *Ii,int n,double rbPar[3])
{
	if(n<=0)return false;
	FittinfPar Par;	
	for(size_t i=0;i<n;i++)	
	{
		Par.epsi.push_back(ei[i]);
		Par.yi.push_back(Ii[i]);	
		//Par.weight.push_bakc(1);
	}	
	rbPar[0]=Fitting(&Par);		
	rbPar[1]=Par.a;
	rbPar[2]=Par.b;	
	/*rbPar[3]=Par.theta;*/
	int iter=0;	
	do
	{
		double tmpL=Fitting_W(&Par);
		double deltL=tmpL-rbPar[0];
		double err=fabs(rbPar[0])+fabs(tmpL);
		if(err>1e-4) err=2*fabs(deltL/err);		
		double pErr=Par.ParError(rbPar[1],rbPar[2]);
		rbPar[0]=tmpL;
		rbPar[1]=Par.a;
		rbPar[2]=Par.b;	
		if(err<0.01||pErr<0.01)	break;		
		iter++;
		//printf("iter #%d\n",iter);
		//printf("LikeHood=%lf,a=%lf,b=%lf\n",rbPar[0],rbPar[1],rbPar[2]);
		//printf("Lerr=%lf,pErr=%lf\n",err,pErr);
	}while(iter<200);
	printf("Iter end, LikeHood=%lf,a=%lf,b=%lf\n",rbPar[0],rbPar[1],rbPar[2]);
	return true;
}



bool IntR_ModelBuilding(vector<double> &ei,vector<double> &Ii,double rbPar[3])
{
	size_t n=ei.size();
	if(n<=0)return false;
	FittinfPar Par;
	for(size_t i=0;i<n;i++)	
	{
		Par.epsi.push_back(ei[i]);
		Par.yi.push_back(Ii[i]);			
	}	
	rbPar[0]=Fitting(&Par);		
	rbPar[1]=Par.a;
	rbPar[2]=Par.b;
	//rbPar[3]=Par.theta;
	int iter=0;	
	do
	{
		double tmpL=Fitting_W(&Par);
		double deltL=tmpL-rbPar[0];
		double err=fabs(rbPar[0])+fabs(tmpL);
		if(err>1e-4) err=2*fabs(deltL/err);
		double pErr=Par.ParError(rbPar[1],rbPar[2]);
		rbPar[0]=tmpL;
		rbPar[1]=Par.a;
		rbPar[2]=Par.b;	
		if(err<0.01||pErr<0.01)	break;
		iter++;
		//printf("iter #%d\n",iter);
		//printf("LikeHood=%lf,a=%lf,b=%lf\n",rbPar[0],rbPar[1],rbPar[2]);
		//printf("Lerr=%lf,pErr==%lf\n",err,pErr);
	}while(iter<100);
	printf("Iter end, LikeHood=%lf,a=%lf,b=%lf\n",rbPar[0],rbPar[1],rbPar[2]);
	return true;
}

double my_mean(double *x,size_t n)
{
	double m=0;
	for(size_t i=0;i<n;i++)	m+=x[i];
	if(n>0) m/=n;
	return m;
}

double my_mean(vector<double> &x)
{
	double m=0;
	size_t n=x.size();
	for(size_t i=0;i<n;i++)	m+=x[i];
	if(n>0) m/=n;
	return m;
}

double my_STD(double *x,double m,size_t n)
{
	double s=0;
	for(size_t i=0;i<n;i++)	s+=(x[i]-m)*(x[i]-m);
	if(n>0) s/=n;
	s=sqrt(s);
	return s;
}

double my_STD(vector<double> &x,double m)
{
	double s=0;
	size_t n=x.size();
	for(size_t i=0;i<n;i++)	s+=(x[i]-m)*(x[i]-m);
	if(n>0) s/=n;
	s=sqrt(s);
	return s;
}

double my_STD(vector<double> &x)
{
	double m=my_mean(x);
	double s=0;
	size_t n=x.size();
	for(size_t i=0;i<n;i++)	s+=(x[i]-m)*(x[i]-m);
	if(n>0) s/=n;
	s=sqrt(s);
	return s;
}


double my_min(double *x,size_t n)
{
	if(n<=0) return 0;
	double min_x=x[0];
	for(size_t i=1;i<n;i++)
	{
		if(min_x>x[i]) min_x=x[i];
	}
	return min_x;
}

double my_min(vector<double> &x)
{
	size_t n=x.size();
	if(n<=0) return 0;
	double min_x=x[0];
	for(size_t i=1;i<n;i++)
	{
		if(min_x>x[i]) min_x=x[i];
	}
	return min_x;
}

double my_max(double *x,size_t n)
{
	if(n<=0) return 0;
	double max_x=x[0];
	for(size_t i=1;i<n;i++)
	{
		if(max_x<x[i]) max_x=x[i];
	}
	return max_x;
}

double my_max(vector<double> &x)
{
	size_t n=x.size();
	if(n<=0) return 0;
	double max_x=x[0];
	for(size_t i=1;i<n;i++)
	{
		if(max_x<x[i]) max_x=x[i];
	}
	return max_x;
}

int my_count(double *x,double min_lim,double max_lim,size_t n)
{
	int C=0;
	for(size_t i=0;i<n;i++)
	{
		if(x[i]>min_lim&&x[i]<max_lim) C++;
	}
	return C;
}


int my_count(vector<double> &x,double min_lim,double max_lim)
{
	int C=0;
	size_t n=x.size();
	for(size_t i=0;i<n;i++)
	{
		if(x[i]>min_lim&&x[i]<max_lim) C++;
	}
	return C;
}


double my_fit_linear(vector<double> &x,vector<double> &y,double c[2])
{
	double *xt,*yt;
	size_t i,num=x.size();
	if(num<=0) return 0;
	xt=new double[num];
	yt=new double[num];
	double my=0;
	for(i=0;i<num;i++) 
	{
		xt[i]=x[i];
		yt[i]=y[i];		
		my+=yt[i];
	}
	my/=num;
	double cov[3];
	double sumsq;
	gsl_fit_linear(xt,1,yt,1,num,c,c+1,cov,cov+1,cov+2,&sumsq);
	double yvs=0;
	for(i=0;i<num;i++)
	{
		yvs+=(yt[i]-my)*(yt[i]-my);		
	}	
	sumsq=1-sumsq/yvs;
	delete []xt;
	delete []yt;
	return sumsq;
}

double getpValue(double x,double mu,double a,double b,double ei)
{
	double sigma=a*exp(-0.5*x)+b;
	ei-=mu;
	double p=0;
	ei=fabs(ei);
	p=gsl_cdf_gaussian_P(ei,sigma);
	p=1-p;	
	return p;
}

double my_sum(double *x,size_t n)
{
	double s=0;
	for(size_t i=0;i<n;i++) s+=x[i];
	return s;
}

double my_sum(vector<double> &x)
{
	double s=0;
	size_t n=x.size();
	for(size_t i=0;i<n;i++) s+=x[i];
	return s;
}


bool my_RobustNormalFit(double *x,size_t n,double *mu,double *sigma,double *PW)
{
	*mu=my_mean(x,n);
	*sigma=my_STD(x,*mu,n);
	double a=my_min(x,n);
	double b=my_max(x,n);
	double alim=*mu-2*(*sigma);
	double blim=*mu+2*(*sigma);
	double L=alim-a+b-blim;
	double K=b-a;
	if(L<=0) return false; 
	size_t deltn=n-my_count(x,alim,blim,n);	
	double PW2=K*deltn/(L*n);
	double PW1=1-PW2;
	double muold=*mu;
	double sigmaold=*sigma;
	double e=1;
	int iter=0;
	double *PWP1,*PWAll;
	PWP1=new double[n];
	PWAll=new double[n];
	size_t i;
	while(e>0.0001)  
	{
		for(i=0;i<n;i++)
		{
			PWP1[i]=PW1*normpdf(x[i],muold,sigmaold);
			PWAll[i]=PWP1[i]+PW2/K;
			PWP1[i]=PWP1[i]/PWAll[i];
		}		
		PW1=my_sum(PWP1,n);
		*mu=0;
		for(i=0;i<n;i++)*mu+=PWP1[i]*x[i];		
		*mu/=PW1;
		*sigma=0;
		for(i=0;i<n;i++)*sigma+=PWP1[i]*(x[i]-muold)*(x[i]-muold);
		*sigma=sqrt(*sigma/PW1);
		PW1=PW1/n;
		PW2=1-PW1;
		e=fabs(*sigma-sigmaold)+fabs(*mu-muold);
		muold=*mu;
		sigmaold=*sigma;
		iter=iter+1;
		if(iter>100) break;
	}
	delete []PWP1;
	delete []PWAll;
	return true;
}

bool my_RobustNormalFit(vector<double> &x,double *mu,double *sigma,double *PW)
{
	size_t n=x.size();
	*mu=my_mean(x);
	*sigma=my_STD(x,*mu);
	double a=my_min(x);
	double b=my_max(x);
	double alim=*mu-2*(*sigma);
	double blim=*mu+2*(*sigma);
	double L=alim-a+b-blim;
	double K=b-a;
	if(L<=0) return false; 
	size_t deltn=n-my_count(x,alim,blim);	
	double PW2=K*deltn/(L*n);
	double PW1=1-PW2;
	double muold=*mu;
	double sigmaold=*sigma;
	double e=1;
	int iter=0;
	double *PWP1,*PWAll;
	PWP1=new double[n];
	PWAll=new double[n];
	size_t i;
	while(e>0.0001)  
	{
		for(i=0;i<n;i++)
		{
			PWP1[i]=PW1*normpdf(x[i],muold,sigmaold);
			PWAll[i]=PWP1[i]+PW2/K;
			PWP1[i]=PWP1[i]/PWAll[i];
		}		
		PW1=my_sum(PWP1,n);
		*mu=0;
		for(i=0;i<n;i++)*mu+=PWP1[i]*x[i];		
		*mu/=PW1;
		*sigma=0;
		for(i=0;i<n;i++)*sigma+=PWP1[i]*(x[i]-muold)*(x[i]-muold);
		*sigma=sqrt(*sigma/PW1);
		PW1=PW1/n;
		PW2=1-PW1;
		e=fabs(*sigma-sigmaold)+fabs(*mu-muold);
		muold=*mu;
		sigmaold=*sigma;
		iter=iter+1;
		if(iter>100) break;
	}
	delete []PWP1;
	delete []PWAll;
	return true;
}
double normpdf(double x,double mu,double sigma)
{
	if(sigma<=0) return 0;
	double t=(x-mu)/sigma;
	t=-t*t;
	return PIS*exp(t/2)/sigma;
}

void normpdfV(double *x,size_t n,double *P,double mu,double sigma)
{
	if(sigma<=0) return;
	for(size_t i=0;i<n;i++)
	{
		double t=(x[i]-mu)/sigma;
		t=-t*t;
		P[i]=PIS*exp(t/2)/sigma;
	}
}

void normpdfV(vector<double> &x,vector<double> &P,double mu,double sigma)
{
	if(sigma<=0) return;
	P.clear();
	size_t n=x.size();
	for(size_t i=0;i<n;i++)
	{
		double t=(x[i]-mu)/sigma;
		t=-t*t;
		P.push_back(PIS*exp(t/2)/sigma);
	}
}

bool IntR_ModelBuildingNew(vector<double> &epsiS,vector<double> &yiS,double rbPar[3])
{
	vector<double> yiW;
	vector<double>sdiW;
	size_t i,sg=yiS.size();
	for(i=0;i<sg;i++) yiS[i]=sqrt(yiS[i]);

	double y_min=my_min(yiS);
	double y_max=my_max(yiS);
	double yRange=y_max-y_min;	
	double yw=yRange/5;
	double ystep=yRange/50;
	double yb=y_min;
	double ye=y_min+yw;
	double yi=y_min+yw/2;
	size_t Ty=yiS.size();
	vector<double> tmpEPS;
	double mu,sigma,pw;
	OutlierRem ot;
	while(yb<y_max)
	{
		for(i=0;i<Ty;i++)
		{
			if(yiS[i]>=yb&&yiS[i]<ye) tmpEPS.push_back(epsiS[i]);
		}
		
		if(tmpEPS.size()<=10) break;
		mu=0;
		sigma=0;
		pw=0;		
		if(ot.oRemoveGTest(tmpEPS))sdiW.push_back(ot.std_new);	
		else sdiW.push_back(my_STD(tmpEPS));
		yiW.push_back(yi);
		tmpEPS.clear();
		yi+=ystep;
		yb+=ystep;
		ye+=ystep;
	}

	Ty=yiW.size();
	if(Ty<=5) 
	{
		printf("modeling failure,check your data.\n");
		return 0;
	}	
	for(i=0;i<Ty;i++) yiW[i]=1/yiW[i];

	double C[2];
	sdiW[0]=1.12*sdiW[0];
		
	rbPar[0]=my_fit_linear(yiW,sdiW,C);	
	rbPar[1]=C[1];
	rbPar[2]=C[0];
	printf("Iter end, LikeHood=%lf,a=%lf,b=%lf\n",rbPar[0],rbPar[1],rbPar[2]);
	return true;
}

bool IntR_ModelBuildingNew(double *ei,double *Ii,int n,double rbPar[3])
{

	vector<double>epsiS;
	vector<double>yiS;
	for(int i=0;i<n;i++)
	{
		epsiS.push_back(ei[i]);
		yiS.push_back(Ii[i]);
	}	
	return IntR_ModelBuildingNew(epsiS,yiS,rbPar);
}

//to debug the alg
bool OutPutModelData(char *fname,vector<double> &ei,vector<double> &Ii)
{
	FILE *fp;
	fp=fopen(fname,"w");
	if(fp==NULL) return false;
	size_t n=ei.size();
	for(size_t i=0;i<n;i++)
	{
		fprintf(fp,"%lf\t%lf\n",ei[i],Ii[i]);
	}
	fclose(fp);
	return true;
}

bool OutPutModelData(char *fname,double *ei,double *Ii,int n)
{
	FILE *fp;
	fp=fopen(fname,"w");
	if(fp==NULL) return false;	
	for(int i=0;i<n;i++)
	{
		fprintf(fp,"%lf\t%lf\n",ei[i],Ii[i]);
	}
	fclose(fp);
	return true;
}
