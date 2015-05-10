// GaussFit.cpp: implementation of the GaussFit class.
//
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "GaussFit.h"
/////////////////////////S_G smoothing//////////////////////////////////////////
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
//using namespace std;

void gsl_diff(double *x,double *diff,int n)
{
	diff[0]=x[0];	
	for(int i=1;i<n;i++)
	{
		diff[i]=x[i]-x[i-1];	
	}	
}

double gsl_pow(double x,int j)
{
	double bReturn=0;
	if(x>0)
	{
		bReturn=exp(j*log(x));
	}
	else if(x<0)
	{
		bReturn=exp(j*log(-x));
		if(j%2!=0) bReturn=-bReturn;	
	}
	return bReturn;	
}
// f windows, k order
bool sgolay(double *x,double *y,int n,int f,int k)
{
	if(n<=0) return false;
	f=f>n?n:f;
	f=f-(f-1)%2;
	double *diffx;
	diffx=new double[n];	
	gsl_diff(x,diffx,n);
	double *c;
	c=new double[n];
	int i;
	for(i=0;i<n;i++) c[i]=y[i];	
	int L,R;
	for(i=0;i<n;i++)
	{
		if(i>0&&(x[i]==x[i-1]))
		{
			c[i]=c[i-1];
			continue;
		}
		//matlab style all 1 change with 0,all n change with n-1 
		L=i;R=i;
		while((R<n-1)&&(x[R+1]==x[i]))R++;
		while(L>0&&(x[L-1]==x[i]))L--;
		double tmp=f-(R-L+1)/2.0;
		int HF=(int)(tmp>0?tmp:0);
		// L = min(n-f+1,max(1,L-HF)); % find leftmost point needed
		int tmpi=L-HF;
		L=n-f;//c style	
		tmpi=tmpi>0?tmpi:0;
		L=L>tmpi?tmpi:L;
		while(L>0&&(x[L-1]==x[L]))L--;
		// R = min(n,max(R+HF,L+f-1)); % find rightmost point needed
		R=R+HF;
		tmpi=L+f-1;//c style?
		tmpi=tmpi>R?tmpi:R;
		R=(n-1)>tmpi?tmpi:(n-1);
		while((R<n-1)&&(x[R]==x[R+1])) R++;
		//R is assured to be larger than L
		int m=R-L+1;
		double *subx,*suby;
		subx=new double[m];
		suby=new double[m];
		//keep i unchange,keep k unchange keep f unchange,keep n unchange
		int j;
		for(j=0;j<m;j++)
		{
			subx[j]=x[j+L]-x[i];
			suby[j]=y[j+L];
		}
		double *dq;
		dq=new double[m];
		gsl_diff(subx,dq,m);
		int vrank=1;
		for(j=1;j<m;j++)
		{
			if(dq[j]>0) vrank++;
		}
		int ncols=(k+1)>vrank?vrank:(k+1);
		gsl_matrix *v;
		if(m==ncols) v=gsl_matrix_alloc(m+1,ncols);
		else v=gsl_matrix_alloc(m,ncols);
		//set first col =1
		for(j=0;j<m;j++)
		{
			gsl_matrix_set(v,j,0,1);
		}
		for(j=1;j<ncols;j++)
		{
			for(int row=0;row<m;row++)
			{
				gsl_matrix_set(v,row,j,gsl_pow(subx[row],j));
			}
		} 
		//% Square v may give infs in the \ solution, so force least squares
		//slove v*c=suby;use multi least squares	
		if(m==ncols) 
		{
			for(j=0;j<ncols;j++)gsl_matrix_set(v,m,j,0);		
		}
		gsl_matrix *cov;
		gsl_vector *fit_y,*d;
		d=gsl_vector_alloc(ncols);
		cov=gsl_matrix_alloc(ncols,ncols);
		if(m==ncols)fit_y=gsl_vector_alloc(m+1);
		else fit_y=gsl_vector_alloc(m);
		for(j=0;j<m;j++)
		{
			gsl_vector_set(fit_y,j,suby[j]);
		}
		if(m==ncols)
		{
			gsl_vector_set(fit_y,m,0);
			m++;
		}	
		double chisq;
		gsl_multifit_linear_workspace * work=gsl_multifit_linear_alloc(m,ncols);
		gsl_multifit_linear(v,fit_y,d,cov,&chisq,work);		
		c[i]=gsl_vector_get(d,0);
		delete []dq;
		delete []subx;
		delete []suby;
		gsl_vector_free(d);
		gsl_matrix_free(v);	
		gsl_vector_free(fit_y);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(work);	
	}
	for(i=0;i<n;i++) y[i]=c[i];
	delete []c;
	delete []diffx;
	return true;
}

int gsl_fit_max(double *y,int n)
{
	double ymax=0;
	int idx=-1;
	for(int i=0;i<n;i++)
	{
		if(ymax<y[i])
		{
			ymax=y[i];
			idx=i;
		}
	}
	return idx;
}

double mgsl_fit_min(double *y,int n)
{
	double ymin=MAX_DOUBLE;
	for(int i=0;i<n;i++)
	{
		if(ymin>y[i])ymin=y[i];
	}
	return ymin;
}

bool Initial_c(FIT_DATA *d,double c[4])
{
	double *y=d->y;
	int n=d->n;
	int idx=gsl_fit_max(y,n);
	if(idx<0||idx>n) return false;
	c[0]=y[idx];
	c[1]=d->x[idx];
	int i;
	for(i=idx;i<n;i++)
	{
		if(y[i]<=0.5*c[0]) break;
	}
	int low,high;
	high=i;
	for(i=idx;i>=0;i--)
	{
		if(y[i]<=0.5*c[0]) break;
	}
	low=i;	
	if(low==-1) low=0;
	if(high==n)high=n-1;
	c[2]=d->x[high]-d->x[low];
	c[2]=c[2]/sqrt(log(16.0));
	c[3]=0.01*c[0];//mgsl_fit_min(y,n);
	return true;
}

double Gauss_F(double x, double *beta)
{
	double t=(x-beta[1])/beta[2];
	double Yi=beta[0]*exp(-t*t)+beta[3];
	return Yi;
}

double gsl_norm(double *v,int n)
{
	double sum=0;
	for(int i=0;i<n;i++)
	{
		sum+=v[i]*v[i];
	}
	return sqrt(sum);
}

double gsl_sum(double *v,int n)
{
	double sum=0;
	for(int i=0;i<n;i++) sum+=v[i];
	return sum;
}

int gsl_divide(gsl_matrix *A,double *b,double *x)
{
	double chisq;
	size_t m=A->size1;
	size_t p=A->size2;
	gsl_vector *y,*d;
	gsl_matrix *cov;
	cov=gsl_matrix_alloc(p,p);
	y=gsl_vector_alloc(m);
	d=gsl_vector_alloc(p);
	size_t i;
	for(i=0;i<m;i++) gsl_vector_set(y,i,b[i]);
	gsl_multifit_linear_workspace * work=gsl_multifit_linear_alloc(m,p);
	gsl_multifit_linear(A,y,d,cov,&chisq,work);	
	for(i=0;i<p;i++)x[i]=gsl_vector_get(d,i);
	gsl_vector_free(d);
	gsl_vector_free(y);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(work);
	return (int)p;
}

double gsl_update_beta(gsl_matrix *J,double lambda,double *rplus,double *beta)
{
	size_t p=J->size2;
	size_t n=J->size1;
	double sum;
	double tmpf;
	gsl_matrix *Jplus;
	Jplus=gsl_matrix_alloc(n+p,p);
	double *step;
	step=new double[p];
	double *tmpv;
	tmpv=new double[p];
	size_t i,k;
	for(i=0;i<p;i++)
	{
		sum=0;	
		for(k=0;k<n;k++)
		{
			tmpf=gsl_matrix_get(J,k,i);
			sum+=tmpf*tmpf;
		}
		tmpv[i]=sqrt(lambda*sum);
	}
	for(i=0;i<p;i++)
	{
		for(k=0;k<n;k++)
		{
			tmpf=gsl_matrix_get(J,k,i);	
			gsl_matrix_set(Jplus,k,i,tmpf);
		}
		for(k=n;k<n+p;k++)
		{
			if(k-n==i) tmpf=tmpv[i];
			else tmpf=0;
			gsl_matrix_set(Jplus,k,i,tmpf);	
		}
	}
	// step = Jplus \ rplus;
	gsl_divide(Jplus,rplus,step);
	// beta(:) = beta(:) + step;
	for(i=0;i<p;i++) beta[i]=beta[i]+step[i];
	sum=(int)gsl_norm(step,p);
	delete []step;
	delete []tmpv;
	gsl_matrix_free(Jplus);
	return sum;
}


int gsl_nlfitOLD(double *x,double *y,size_t n,double *beta,size_t p,F f,OP *options)
{
	if(n<=0) return R_ERROR;
	if(p<=0) return R_ERROR;
	double *yfit,*r;
	gsl_matrix *J1;
	J1=NULL;
	yfit=new double[n];
	r=new double[n];
	size_t i;
	for( i=0;i<n;i++)yfit[i]=f(x[i],beta);
	double sqrteps=sqrt(DEPS);	
	J1=gsl_matrix_alloc(n,p);
	gsl_matrix_set_zero(J1);
	for(i=0;i<n;i++) r[i]=y[i]-yfit[i];
	double sse=0;
	for(i=0;i<n;i++) 
	{
		if(y[i]>1e-6||y[i]<-1e-6)
			sse+=r[i]*r[i]/y[i];
	}
	int maxiter=options->MaxIter;
	double betatol = options->TolX;
	double rtol = options->TolFun;
	double fdiffstep = options->DerivStep;

	// Set initial weight for LM algorithm.
	double lambda=0.01;

	int iter=0;
	bool breakOut=false;
	double *betaold;
	betaold=new double[p];
	double sseold;
	double *tmpv;
	tmpv=new double[p];
	double *dy;
	dy=new double[n];
	size_t k;
	double *rplus;
	rplus=new double[n+p];
	double norm_step;
	int bReturn;
	double delta;

	while(iter<maxiter)
	{
		iter++;
		for(i=0;i<p;i++)betaold[i]=beta[i];
		sseold = sse;
		//  Compute a finite difference approximation to the Jacobian
		for(k=0;k<p;k++)
		{
			if (beta[k]==0)
			{
				double nb=sqrt(gsl_norm(beta,p));
				if(nb<DEPS) nb=1;
				delta=fdiffstep*nb;
			}
			else delta=fdiffstep*beta[k];
			for(i=0;i<p;i++)tmpv[i]=beta[i];
			tmpv[k]+=delta;
			for(i=0;i<n;i++)
			{			
				dy[i]=f(x[i],tmpv)-yfit[i];
				gsl_matrix_set(J1,i,k,dy[i]/delta);
			}
		}
		// Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r		
		//Jplus = [J; diag(sqrt(lambda*sum(J.^2)))];
		//rplus = [r; zerosp];
		for(i=0;i<n;i++)rplus[i]=r[i];
		for(i=n;i<n+p;i++)rplus[i]=0;
		norm_step=gsl_update_beta(J1,lambda,rplus,beta);
		// Evaluate the fitted values at the new coefficients and
		// compute the residuals and the SSE.
		//yfit = model(beta,X);
		for(i=0;i<n;i++)yfit[i]=f(x[i],beta);
		//r = y(:) - yfit(:);
		for(i=0;i<n;i++)r[i]=y[i]-yfit[i];
		//sse = r'*r;
		sse=0;
		for(i=0;i<n;i++) 
		{
			if(y[i]>1e-6||y[i]<-1e-6)
				sse+=r[i]*r[i]/y[i];
		}
		//If the LM step decreased the SSE, decrease lambda to downweight the
		// steepest descent direction.
		if(sse<sseold)lambda=0.1*lambda;		
		//If the LM step increased the SSE, repeatedly increase lambda to
		// upweight the steepest descent direction and decrease the step size
		// until we get a step that does decrease SSE.
		else
		{
			 while(sse>sseold)
			 {
				 lambda=10*lambda; 
				 //warning('stats:nlinfit:UnableToDecreaseSSE', ...
				 //'Unable to find a step that will decrease SSE.  Returning results from last iteration.');
				 if(lambda>1e16)
				 {
					 breakOut=true;
					 break;
				 }
				 //Jplus = [J; diag(sqrt(lambda*sum(J.^2)))];
				 //step = Jplus \ rplus; 
				 //beta(:) = betaold(:) + step;
				 for(i=0;i<p;i++) beta[i]=betaold[i];
				 norm_step=gsl_update_beta(J1,lambda,rplus,beta);	
				 //yfit = model(beta,X);
				 for(i=0;i<n;i++)yfit[i]=f(x[i],beta);	
				 //r = y(:) - yfit(:);
				 for(i=0;i<n;i++)r[i]=y[i]-yfit[i];	
				 //sse = r'*r;
				 sse=0;
				 for(i=0;i<n;i++)
				{
					if(y[i]>1e-6||y[i]<-1e-6)
						sse+=r[i]*r[i]/y[i];
				 }
			 }
		}

		//if(options->Disp==2)
	///	{
		//	cout<<iter<<"\t"<<sse<<endl;
		//}

		if(norm_step< betatol*(sqrteps+gsl_norm(beta,p)))
		{
			bReturn=R_MIN_TolX;
			break;	
		}
        //disp('Iterations terminated: relative norm of the current step 
		//is less than OPTIONS.TolX');
		else if(fabs(sse-sseold)<=rtol*sse)
		{
			bReturn=R_MIN_TolFun;
			break;
		}
       //disp('Iterations terminated: relative change in SSE less than OPTIONS.TolFun');
		else if(breakOut)
		{
			bReturn=R_INC_SSE;
			break; 
		}
	}
	if(iter>=maxiter)bReturn=R_MAX_Iter;
		//warning('stats:nlinfit:IterationLimitExceeded', ...
		//'Iteration limit exceeded.  Returning results from final iteration.');
	double sst=0;
	for(i=0;i<n;i++)
	{
		delta+=y[i];
	}
	delta/=n;
	for(i=0;i<n;i++)
	{
		sst+=(y[i]-delta)*(y[i]-delta);
	}
	options->sse=sse;
	options->r=1-sse/sst;
	delete []yfit;
	delete []r;
	delete []rplus;
	delete []betaold;
	delete []tmpv;
	delete []dy;
	gsl_matrix_free(J1);
	return bReturn;
}

int gsl_nlfit(double *x,double *y,size_t n,double *beta,size_t p,F f,OP *options)
{
	if(n<=0) return R_ERROR;
	if(p<=0) return R_ERROR;
	double *yfit,*r;
	gsl_matrix *J1;
	J1=NULL;
	yfit=new double[n];
	r=new double[n];
	size_t i;
	for( i=0;i<n;i++)yfit[i]=f(x[i],beta);
	double sqrteps=sqrt(DEPS);	
	J1=gsl_matrix_alloc(n,p);
	gsl_matrix_set_zero(J1);
	for(i=0;i<n;i++) r[i]=y[i]-yfit[i];
	double sse=0;
	for(i=0;i<n;i++) sse+=r[i]*r[i];
	int maxiter=options->MaxIter;
	double betatol = options->TolX;
	double rtol = options->TolFun;
	double fdiffstep = options->DerivStep;

	// Set initial weight for LM algorithm.
	double lambda=0.01;

	int iter=0;
	bool breakOut=false;
	double *betaold;
	betaold=new double[p];
	double sseold;
	double *tmpv;
	tmpv=new double[p];
	double *dy;
	dy=new double[n];
	size_t k;
	double *rplus;
	rplus=new double[n+p];
	double norm_step;
	int bReturn;
	double delta;

	while(iter<maxiter)
	{
		iter++;
		for(i=0;i<p;i++)betaold[i]=beta[i];
		sseold = sse;
		//  Compute a finite difference approximation to the Jacobian
		for(k=0;k<p;k++)
		{
			if (beta[k]==0)
			{
				double nb=sqrt(gsl_norm(beta,p));
				if(nb<DEPS) nb=1;
				delta=fdiffstep*nb;
			}
			else delta=fdiffstep*beta[k];
			for(i=0;i<p;i++)tmpv[i]=beta[i];
			tmpv[k]+=delta;
			for(i=0;i<n;i++)
			{			
				dy[i]=f(x[i],tmpv)-yfit[i];
				gsl_matrix_set(J1,i,k,dy[i]/delta);
			}
		}
		// Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r		
		//Jplus = [J; diag(sqrt(lambda*sum(J.^2)))];
		//rplus = [r; zerosp];
		for(i=0;i<n;i++)rplus[i]=r[i];
		for(i=n;i<n+p;i++)rplus[i]=0;
		norm_step=gsl_update_beta(J1,lambda,rplus,beta);
		// Evaluate the fitted values at the new coefficients and
		// compute the residuals and the SSE.
		//yfit = model(beta,X);
		for(i=0;i<n;i++)yfit[i]=f(x[i],beta);
		//r = y(:) - yfit(:);
		for(i=0;i<n;i++)r[i]=y[i]-yfit[i];
		//sse = r'*r;
		sse=0;
		for(i=0;i<n;i++) sse+=r[i]*r[i];
		//If the LM step decreased the SSE, decrease lambda to downweight the
		// steepest descent direction.
		if(sse<sseold)lambda=0.1*lambda;		
		//If the LM step increased the SSE, repeatedly increase lambda to
		// upweight the steepest descent direction and decrease the step size
		// until we get a step that does decrease SSE.
		else
		{
			 while(sse>sseold)
			 {
				 lambda=10*lambda; 
				 //warning('stats:nlinfit:UnableToDecreaseSSE', ...
				 //'Unable to find a step that will decrease SSE.  Returning results from last iteration.');
				 if(lambda>1e5)
				 {
					 breakOut=true;
					 break;
				 }
				 //Jplus = [J; diag(sqrt(lambda*sum(J.^2)))];
				 //step = Jplus \ rplus; 
				 //beta(:) = betaold(:) + step;
				 for(i=0;i<p;i++) beta[i]=betaold[i];
				 norm_step=gsl_update_beta(J1,lambda,rplus,beta);	
				 //yfit = model(beta,X);
				 for(i=0;i<n;i++)yfit[i]=f(x[i],beta);	
				 //r = y(:) - yfit(:);
				 for(i=0;i<n;i++)r[i]=y[i]-yfit[i];	
				 //sse = r'*r;
				 sse=0;
				 for(i=0;i<n;i++) sse+=r[i]*r[i];
			 }
		}

		//if(options->Disp==2)
	///	{
		//	cout<<iter<<"\t"<<sse<<endl;
		//}

		if(norm_step< betatol*(sqrteps+gsl_norm(beta,p)))
		{
			bReturn=R_MIN_TolX;
			break;	
		}
        //disp('Iterations terminated: relative norm of the current step 
		//is less than OPTIONS.TolX');
		else if(fabs(sse-sseold)<=rtol*sse)
		{
			bReturn=R_MIN_TolFun;
			break;
		}
       //disp('Iterations terminated: relative change in SSE less than OPTIONS.TolFun');
		else if(breakOut)
		{
			bReturn=R_INC_SSE;
			break; 
		}
	}
	if(iter>=maxiter)bReturn=R_MAX_Iter;
		//warning('stats:nlinfit:IterationLimitExceeded', ...
		//'Iteration limit exceeded.  Returning results from final iteration.');
	double sst=0;
	for(i=0;i<n;i++)
	{
		delta+=y[i];
	}
	delta/=n;
	for(i=0;i<n;i++)
	{
		sst+=(y[i]-delta)*(y[i]-delta);
	}
	options->sse=sse;
	options->r=1-sse/sst;
	delete []yfit;
	delete []r;
	delete []rplus;
	delete []betaold;
	delete []tmpv;
	delete []dy;
	gsl_matrix_free(J1);
	return bReturn;
}


int gsl_nlfitNew(double *x,double *y,size_t n,double *beta,size_t p,FNew f,OP *options)
{
	if(n<=0) return R_ERROR;
	if(p<=0) return R_ERROR;
	if(p%4!=0) return R_ERROR;
	int num=p/4;
	double *yfit,*r;
	gsl_matrix *J1;
	J1=NULL;
	yfit=new double[n];
	r=new double[n];
	size_t i;
	for( i=0;i<n;i++)yfit[i]=f(x[i],beta,num);
	double sqrteps=sqrt(DEPS);	
	J1=gsl_matrix_alloc(n,p);
	gsl_matrix_set_zero(J1);
	for(i=0;i<n;i++) r[i]=y[i]-yfit[i];
	double sse=0;
	for(i=0;i<n;i++) sse+=r[i]*r[i];
	int maxiter=options->MaxIter;
	double betatol = options->TolX;
	double rtol = options->TolFun;
	double fdiffstep = options->DerivStep;

	// Set initial weight for LM algorithm.
	double lambda=0.01;

	int iter=0;
	bool breakOut=false;
	double *betaold;
	betaold=new double[p];
	double sseold;
	double *tmpv;
	tmpv=new double[p];
	double *dy;
	dy=new double[n];
	size_t k;
	double *rplus;
	rplus=new double[n+p];
	double norm_step;
	int bReturn;
	double delta;

	while(iter<maxiter)
	{
		iter++;
		for(i=0;i<p;i++)betaold[i]=beta[i];
		sseold = sse;
		//  Compute a finite difference approximation to the Jacobian
		for(k=0;k<p;k++)
		{
			if (beta[k]==0)
			{
				double nb=sqrt(gsl_norm(beta,p));
				if(nb<DEPS) nb=1;
				delta=fdiffstep*nb;
			}
			else delta=fdiffstep*beta[k];
			for(i=0;i<p;i++)tmpv[i]=beta[i];
			tmpv[k]+=delta;
			for(i=0;i<n;i++)
			{			
				dy[i]=f(x[i],tmpv,num)-yfit[i];
				gsl_matrix_set(J1,i,k,dy[i]/delta);
			}
		}
		// Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r		
		//Jplus = [J; diag(sqrt(lambda*sum(J.^2)))];
		//rplus = [r; zerosp];
		for(i=0;i<n;i++)rplus[i]=r[i];
		for(i=n;i<n+p;i++)rplus[i]=0;
		norm_step=gsl_update_beta(J1,lambda,rplus,beta);
		// Evaluate the fitted values at the new coefficients and
		// compute the residuals and the SSE.
		//yfit = model(beta,X);
		for(i=0;i<n;i++)yfit[i]=f(x[i],beta,num);
		//r = y(:) - yfit(:);
		for(i=0;i<n;i++)r[i]=y[i]-yfit[i];
		//sse = r'*r;
		sse=0;
		for(i=0;i<n;i++) sse+=r[i]*r[i];
		//If the LM step decreased the SSE, decrease lambda to downweight the
		// steepest descent direction.
		if(sse<sseold)lambda=0.5*lambda;		
		//If the LM step increased the SSE, repeatedly increase lambda to
		// upweight the steepest descent direction and decrease the step size
		// until we get a step that does decrease SSE.
		else
		{
			 while(sse>sseold)
			 {
				 lambda=2*lambda; 
				 //warning('stats:nlinfit:UnableToDecreaseSSE', ...
				 //'Unable to find a step that will decrease SSE.  Returning results from last iteration.');
				 if(lambda>1e5)
				 {
					 breakOut=true;
					 break;
				 }
				 //Jplus = [J; diag(sqrt(lambda*sum(J.^2)))];
				 //step = Jplus \ rplus; 
				 //beta(:) = betaold(:) + step;
				 for(i=0;i<p;i++) beta[i]=betaold[i];
				 norm_step=gsl_update_beta(J1,lambda,rplus,beta);	
				 //yfit = model(beta,X);
				 for(i=0;i<n;i++)yfit[i]=f(x[i],beta,num);	
				 //r = y(:) - yfit(:);
				 for(i=0;i<n;i++)r[i]=y[i]-yfit[i];	
				 //sse = r'*r;
				 sse=0;
				 for(i=0;i<n;i++) sse+=r[i]*r[i];
			 }
		}

		//if(options->Disp==2)
	///	{
		//	cout<<iter<<"\t"<<sse<<endl;
		//}

		if(norm_step< betatol*(sqrteps+gsl_norm(beta,p)))
		{
			bReturn=R_MIN_TolX;
			break;	
		}
        //disp('Iterations terminated: relative norm of the current step 
		//is less than OPTIONS.TolX');
		else if(fabs(sse-sseold)<=rtol*sse)
		{
			bReturn=R_MIN_TolFun;
			break;
		}
       //disp('Iterations terminated: relative change in SSE less than OPTIONS.TolFun');
		else if(breakOut)
		{
			bReturn=R_INC_SSE;
			break; 
		}
	}
	if(iter>=maxiter)bReturn=R_MAX_Iter;
		//warning('stats:nlinfit:IterationLimitExceeded', ...
		//'Iteration limit exceeded.  Returning results from final iteration.');
	double sst=0;
	for(i=0;i<n;i++)
	{
		delta+=y[i];
	}
	delta/=n;
	for(i=0;i<n;i++)
	{
		sst+=(y[i]-delta)*(y[i]-delta);
	}
	options->sse=sse;
	options->r=1-sse/sst;
	delete []yfit;
	delete []r;
	delete []rplus;
	delete []betaold;
	delete []tmpv;
	delete []dy;
	gsl_matrix_free(J1);
	return bReturn;
}

double Gauss_Fit(FIT_DATA *d,double beta[4])
{
	if(!Initial_c(d,beta)) return 0;
	//gsl_get_sigma(d);
	OP op;
	op.DerivStep=exp(log(DEPS)/3);
	op.Disp=2;
	op.MaxIter=100;
	op.TolFun=1e-8;
	op.TolX=1e-8;
	gsl_nlfit(d->x,d->y,d->n,beta,4,&Gauss_F,&op);
	return op.r;
}

double gsl_gauss(double x, void *params)
{
	double *P;
	P=(double *)params;
	double A=    P[0];
	double mu=   P[1];
	double deta= P[2];
	//double b=    P[3];
	double t=(x-mu)/deta;
	//return A*exp(-t*t)+b;
	return A*exp(-t*t);
}

double GetArea(double beta[3],double a,double b)
{
	gsl_integration_workspace *w=gsl_integration_workspace_alloc(200);
	gsl_function F;
	F.function = &gsl_gauss;
	F.params = (void *)beta;
	double result,error;
	gsl_integration_qags(&F,a,b,0,1e-7,200,w, &result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

bool gsl_GetPShape(double *x,double *y,int n,double *PH,double *PW,double *PA)
{
	if(n<=1) return false;
	gsl_interp_accel *acc=gsl_interp_accel_alloc();
	gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline,n);
	gsl_spline_init(spline,x,y,n);
	double w=(x[n-1]-x[0])/100;
	double xi=x[0];
	double *y1;
	y1=new double[100];
	int i;
	for(i=0;i<100;i++)
	{		
		y1[i]=gsl_spline_eval(spline,xi,acc);
		xi+=w;
	}
	double maxh=0;
	int idx=-1;
	for(i=0;i<100;i++)
	{
		if(maxh<y1[i])
		{
			maxh=y1[i];
			idx=i;
		}
	}
	if(idx==-1)
	{
		delete []y1;
		gsl_spline_free (spline);
		gsl_interp_accel_free(acc);
		return false;
	}
	*PH=maxh;
	for(i=idx;i<100;i++)
	{
		if(y1[i]<maxh*0.5)break;
	}
	if(i==100)i--;
	int h=i;
	for(i=idx;i>=0;i--)
	{
		if(y1[i]<0.5*maxh) break;
	}
	if(i<0) i=0;
	int l=i;
	*PW=w*(h-l);
	*PA=0;
	for(i=1;i<100;i++) *PA+=(y1[i]+y1[i-1])*w/2;
	delete []y1;
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return true;
}

double gsl_CShape(double x,double beta[5])
{
	double H=beta[0];
	double M=beta[1];
	double a=beta[2];
	double D=beta[3];
	double b=beta[4];
	if(a>=2||a<=0) return 0;
	double cut=M-D*(4-a*a)/(2*a);
    if(x<cut) return b;
    else
	{
        double t=2*a*(x-M)/(D*(4-a*a));
        return H*exp((4/(a*a)-1)*(log(1+t)-t))+b;
	}
}

double gsl_CShapeNew(double x,double beta[4])
{
	double H=beta[0];
	double M=beta[1];
	double a=beta[2];
	double D=beta[3];
	if(a>=2||a<=0) return 0;
	double cut=M-D*(4-a*a)/(2*a);
    if(x<cut) return 0;
    else
	{
        double t=2*a*(x-M)/(D*(4-a*a));
        return H*exp((4/(a*a)-1)*(log(1+t)-t));
	}
}


double gsl_ChroFit(double *x,double *y,int n,double beta[5])
{
	if(n<=1) return 0;
	int idx=gsl_fit_max(y,n);
	if(idx<0||idx>n) return 0;
	beta[0]=y[idx];
	beta[1]=x[idx];
	beta[2]=1;
	beta[3]=0.3;
	beta[4]=0.01*beta[0];
	OP op;
	op.DerivStep=exp(log(DEPS)/3);
	op.Disp=2;
	op.MaxIter=100;
	op.TolFun=1e-8;
	op.TolX=1e-8;
	if(gsl_nlfit(x,y,n,beta,5,&gsl_CShape,&op)==R_INC_SSE)
	{
		 return 0;
	}
	return op.r;
}

double w_mean(double *Res, char *ResW, int n)
{
   int num=0;
   double m=0;
   for(int i=0;i<n;i++)
   {
	   if(ResW[i])
	   {
		   m+=Res[i];
		   num++;
	   }
   }
   if(num>0) m/=num;
   return m;
}

double gsl_w_mean(gsl_vector *x,gsl_vector *w)
{
   double m=0;
   int n=x->size;
   for(int i=0;i<n;i++)
   { 
	   m+=gsl_vector_get(x,i)*gsl_vector_get(w,i);	  
   }
   if(n>0) m/=n;
   return m;
}

double gsl_w_std(gsl_vector *x,gsl_vector *w,double m)
{
   double sd=0;
   int n=x->size;
   for(int i=0;i<n;i++)
   { 
	   double xi=gsl_vector_get(x,i);
	   double wi=gsl_vector_get(w,i);
	   xi-=m;
	   sd+=wi*xi*xi;	   	  
   }
   if(n>0) sd/=n;
   return sd;
}

double gsl_mean(gsl_vector *x)
{
   double m=0;
   int n=x->size;
   for(int i=0;i<n;i++)
   { 
	   m+=gsl_vector_get(x,i);	  
   }
   if(n>0) m/=n;
   return m;
}

double w_mean(gsl_vector *Res, char *ResW, int n)
{
   int num=0;
   double m=0;
   for(int i=0;i<n;i++)
   {
	   if(ResW[i])
	   {
		   m+=gsl_vector_get(Res,i);
		   num++;
	   }
   }
   if(num>0) m/=num;
   return m;
}

double w_std(double *Res, char *ResW, int n,double m)
{
	double std=0;
	int num=0;
	for(int i=0;i<n;i++)
	{
	   if(ResW[i])
	   {
		   std+=(Res[i]-m)*(Res[i]-m);
		   num++;
	   }
	}
	if(num>0) std=sqrt(std/num);
	return std;
}

double w_std(gsl_vector *Res, char *ResW, int n,double m)
{
	double std=0;
	int num=0;
	double el;
	for(int i=0;i<n;i++)
	{
	   if(ResW[i])
	   {
		   el=gsl_vector_get(Res,i);
		   std+=(el-m)*(el-m);
		   num++;
	   }
	}
	if(num>0) std=sqrt(std/num);
	return std;
}
///////////functions for robust fitting
int gsl_vector_output(gsl_vector *X)
{
	for(UINT i=0;i<X->size;i++)
	{
		printf("%f\n",gsl_vector_get(X,i));
	}
	return X->size;
}

int gsl_matrix_output(gsl_matrix *X)
{
	int row=X->size1;
	int col=X->size2;
	for(int i=0;i<row;i++)
	{
		for(int k=0;k<col;k++)
		{
			printf("%f\t",gsl_matrix_get(X,i,k));
		}
		printf("\n");	
	}
	return col*row;
}


int gsl_matrix_inverse(gsl_matrix *X)
{
	int row=X->size1;
	int col=X->size2;
	if(row!=col) return 0;
	gsl_permutation *p;
	gsl_matrix *LU;	
	LU=gsl_matrix_alloc(row,col);
	gsl_matrix_memcpy(LU,X);
	p=gsl_permutation_alloc(col);
	int s;
	gsl_linalg_LU_decomp (LU,p,&s);
	gsl_linalg_LU_invert(LU,p,X);
	gsl_permutation_free(p);
	gsl_matrix_free(LU);
	return col;
}

double gsl_vector_sum(gsl_vector *X)
{
	double sum=0;
	for(UINT i=0;i<X->size;i++)
	{
		sum+=gsl_vector_get(X,i);
	}
	return sum;

}
///Y is a square matrix,used by div
int gsl_matrix_mul(gsl_matrix *X,gsl_matrix *Y)
{

	int row=X->size1;
	int col=X->size2;
	int row1=Y->size1;
	int col1=Y->size2;
	int i,k;
	double el;
	if(col!=row1) return 0;
	gsl_vector *x_row,*y_col;
	x_row=gsl_vector_alloc(col);
	y_col=gsl_vector_alloc(col);
	for(i=0;i<row;i++)
	{
		gsl_matrix_get_row(x_row,X,i);
		for(k=0;k<col;k++)
		{
			gsl_matrix_get_col(y_col,Y,k);
			gsl_vector_mul(y_col,x_row);
			el=gsl_vector_sum(y_col);
			gsl_matrix_set(X,i,k,el);
		}
	}
	gsl_vector_free(x_row);
	gsl_vector_free(y_col);
	return col*row;
}

int gsl_matrix_div(gsl_matrix *X,gsl_matrix *Y)
{
	gsl_matrix_inverse(Y);
	int bReturn=gsl_matrix_mul(X,Y);
	return bReturn;	
}

bool gsl_adfactor_get(gsl_matrix *X,gsl_vector *h)
{
	int i,k;
	gsl_matrix *X1,*R;	
	gsl_vector *tau,*tmp;
	int row=X->size1;
	int col=X->size2;
	//the number of observations must be great than the the variables
	if(row<col) return false;
	X1=gsl_matrix_alloc(row,col);
	gsl_matrix_memcpy(X1,X);
	int size_tau=row>col?col:row;
	tau=gsl_vector_alloc(size_tau);
	gsl_linalg_QR_decomp(X1,tau);
	//////////Get the compact format R and X
	R=gsl_matrix_alloc(col,col);	
	tmp=gsl_vector_alloc(col);	
	for(k=0;k<col;k++)
	{
		gsl_matrix_get_row(tmp,X1,k);
		gsl_matrix_set_row(R,k,tmp);
	}	
	for(i=0;i<col;i++)
	{
		for(k=0;k<col;k++)
		{
			if(i>k)	gsl_matrix_set(R,i,k,0);
		}
	}
	gsl_matrix_memcpy(X1,X);
	///////rearrange R and X1,where diag(R) is destend
	double tol=EPS*row;
	double el;
	for(i=0;i<col;i++)
	{	
		el=gsl_matrix_get(R,i,i);
		if(fabs(el)<tol)
		{
			gsl_matrix_set(R,i,i,tol);
		}
	}
	/////////Get X1/R;
	gsl_matrix_div(X1,R);
	//////////Get X1^2
	gsl_matrix_mul_elements(X1,X1);
	for(i=0;i<row;i++)
	{
		gsl_matrix_get_row(tmp,X1,i);
		el=gsl_vector_sum(tmp);
		if(el>0.9999) el=0.9999;
		el=1/sqrt(1-el);
		gsl_vector_set(h,i,el);
	}
	gsl_matrix_free(X1);
	gsl_matrix_free(R);
	gsl_vector_free(tau);
	gsl_vector_free(tmp);
	return true;
}

int gsl_sort_double(double *data,int n)
{
	double el;
	for(int i=0;i<n;i++)
	{
		for(int k=i+1;k<n;k++)
		{
			if(data[i]>data[k])
			{
				el=data[i];
				data[i]=data[k];
				data[k]=el;
			}
		}
	}
	return n;
}

double gsl_stats_median(double *data,int n)
{
	int idx,idx1;
	if(n%2==1)
	{
		idx=(n-1)/2;
		return data[idx];
	}
	else
	{
		idx=n/2;
		idx1=idx-1;
		if(idx1<0) idx1=0;
		return (data[idx]+data[idx1])/2;
	}
}

double gsl_madsigma(gsl_vector *Res,int p)
{
	int n=Res->size;
	if(n<=0) return 0;
	double *e,*rs;
	e=new double[n];
	rs=new double[n];
	int i;
	for(i=0;i<n;i++)
	{
		e[i]=gsl_vector_get(Res,i);
	}	
	gsl_sort_double(e,n);
	double m=gsl_stats_median(e,n);
	for(i=0;i<n;i++)
	{
		rs[i]=fabs(e[i]-m);
	}
	gsl_sort_double(rs,n);
	if(fabs(m)>rs[n-1])
	{
		for(i=0;i<n;i++)
		{
			rs[i]=fabs(e[i]);
		}
	}
	double s=gsl_stats_median(rs+p,n-p);
	s=s/0.6745;
	if(s==0)s=0.5*gsl_stats_mean(rs,1,n);
	delete []rs;
	delete []e;
	return s;
}

int gsl_update_w(gsl_vector *w,gsl_vector *Res,gsl_vector *h,int type,int p)
{
	UINT n=Res->size;
	if(n!=h->size) return 0;
	double el;
	UINT i;
	for(i=0;i<n;i++)
	{
		el=gsl_vector_get(Res,i)*gsl_vector_get(h,i);
		gsl_vector_set(Res,i,el);
	}
	double s=gsl_madsigma(Res,p);
	if(s<=0) s=1;
	double tune;
	double S_EPS=sqrt(EPS);
	switch(type)
	{
	case 0:// case 'andrews',  tune = 1.339; 
		{
			tune=1.339;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);
				if(el<S_EPS)el=S_EPS;
				if(el>PI)el=0;
				else el=sin(el)/el;
				gsl_vector_set(w,i,el);
			}
			break;
		}
	case 1://   case 'bisquare', tune = 4.685;
		{
			tune=4.685;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);			
				if(el>=1)el=0;
				else 
				{
					el=1-el*el;
					el=el*el;
				}
				gsl_vector_set(w,i,el);
			}
			break;			
		}
	case 2://    case 'cauchy',   tune = 2.385;
		{
			tune=2.385;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);
				el=1/(1+el*el);
				gsl_vector_set(w,i,el);
			}
			break;
		}
	case 3://    case 'fair',     tune = 1.400;
		{
			tune=1.400;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);
				el=1/(1+fabs(el));
				gsl_vector_set(w,i,el);
			}
			break;
		}
	case 4://    case 'huber',    tune = 1.345;
		{
			tune=1.345;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);
				if(el<1)el=1;			
				gsl_vector_set(w,i,el);
			}
			break;
		}
	case 5://    case 'logistic', tune = 1.205;
		{
			tune=1.205;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);
				if(el<S_EPS) el=S_EPS;
				el=tanh(el)/el;				
				gsl_vector_set(w,i,el);
			}
			break;
		}
	case 6://    case 'ols',      tune = 1;
		{	
			for(i=0;i<n;i++)
			{			
				gsl_vector_set(w,i,1);
			}
			break;
		}
	case 7://    case 'talwar',   tune = 2.795; 
		{
			tune=2.795;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);
				if(el<1) el=1;
				else el=0;							
				gsl_vector_set(w,i,el);
			}
			break;
		}
	case 8://   case 'welsch',   tune = 2.985;
		{
			tune=2.985;
			for(i=0;i<n;i++)
			{
				el=gsl_vector_get(Res,i);
				el/=(s*tune);
				el=exp(-el*el);						
				gsl_vector_set(w,i,el);
			}
			break;
		}
	default:// same as case 6
		{
			for(i=0;i<n;i++)
			{			
				gsl_vector_set(w,i,1);
			}	
		}
	}
	return n;
}

int gsl_res_get(gsl_matrix *X,gsl_vector *y,gsl_vector *c,gsl_vector *Res)
{
	double el;
	int n=X->size1;
	int p=X->size2;
	for(int i=0;i<n;i++)
	{	
		el=0;
		for(int j=0;j<p;j++)
		{
			el+=gsl_matrix_get(X,i,j)*gsl_vector_get(c,j);
		}
		el=gsl_vector_get(y,i)-el;
		gsl_vector_set(Res,i,el);		
	}
	return n;
}

int gsl_res_get(gsl_matrix *X,gsl_vector *y,gsl_vector *c,double *Res)
{
	double el;
	int n=X->size1;
	int p=X->size2;
	for(int i=0;i<n;i++)
	{	
		el=0;
		for(int j=0;j<p;j++)
		{
			el+=gsl_matrix_get(X,i,j)*gsl_vector_get(c,j);
		}
		Res[i]=gsl_vector_get(y,i)-el;				
	}
	return n;
}

int gsl_count(double *x,double lmin,double lmax,int n)
{
	int c=0;
	for(int i=0;i<n;i++)
	{
		if(x[i]>lmin&&x[i]<lmax) c++;
	}
	return c;
}

#define PIS 0.398942280401433
double gsl_normpdf(double x,double mu,double sigma)
{
	if(sigma<=0) return 0;
	double t=(x-mu)/sigma;
	t=-t*t;
	return PIS*exp(t/2)/sigma;
}

bool gsl_RobustNormalFit(double *x,int n,double *mu,double *sigma,double *PW)
{	
	*mu=gsl_stats_mean(x,1,n);
	*sigma=gsl_stats_sd_m(x,1,n,*mu);
	double a=gsl_stats_min(x,1,n);
	double b=gsl_stats_max(x,1,n);
	double alim=*mu-2*(*sigma);
	double blim=*mu+2*(*sigma);
	double L=alim-a+b-blim;
	double K=b-a;
	if(L<=0) return false; 
	size_t deltn=n-gsl_count(x,alim,blim,n);	
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
			PWP1[i]=PW1*gsl_normpdf(x[i],muold,sigmaold);
			PWAll[i]=PWP1[i]+PW2/K;
			PWP1[i]=PWP1[i]/PWAll[i];
		}		
		PW1=gsl_sum(PWP1,n);
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

bool gsl_RobustNormalFit(gsl_vector *x,double *mu,double *sigma,double *PW)
{	
	int n=x->size;
	if(n<=0) return false;
	double *tx;
	tx=new double[n];
	for(int i=0;i<n;i++) tx[i]=gsl_vector_get(x,i);
	bool bR=gsl_RobustNormalFit(tx,n,mu,sigma,PW);
	delete []tx;
	return bR;
}

double gsl_RSquare(gsl_vector *y,gsl_vector *Res,gsl_vector *w)
{
	int n=y->size;
	double sum1=0;
	double sum2=0;
	double ym=gsl_w_mean(y,w);
	for(int i=0;i<n;i++)
	{
		double r=gsl_vector_get(Res,i);
		double wt=gsl_vector_get(w,i);
		double yt=gsl_vector_get(y,i);
		sum1+=wt*r*r;
		yt-=ym;
		sum2+=yt*yt*wt;
	}	
	if(sum2>0) sum1/=sum2;
	else sum1=0;
	return 1-sum1;
}

double gsl_two_step_reg(double *xi,double *yi,int p,int n,double par[])
{
	gsl_matrix *X;
	gsl_vector *y,*c;	
	if(p<=1) return 0;
	if(n<=p) return 0;
	X=gsl_matrix_alloc(n,p);
	y=gsl_vector_alloc(n);
	c=gsl_vector_alloc(p);
	size_t i,k;
	for(i=0;i<n;i++)
	{
		gsl_vector_set(y,i,yi[i]);
		double tx=1;
		for(k=0;k<p;k++)
		{
			gsl_matrix_set(X,i,k,tx);
			tx*=xi[i];
		}
	}	

	gsl_matrix *cov;
	gsl_vector *w,*Res;
	double chisq;
	w=gsl_vector_alloc(n);
	Res=gsl_vector_alloc(n);
	cov=gsl_matrix_alloc(p,p);
	
	for(i=0;i<n;i++)
	{
		gsl_vector_set(w,i,1.0);
	}
	
	gsl_multifit_linear_workspace * work= gsl_multifit_linear_alloc(n,p);	
	gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);	
	gsl_res_get(X,y,c,Res);
	double mu,sigma,PW;
	gsl_RobustNormalFit(Res,&mu,&sigma,&PW);

	double LL=mu-6*sigma;
	double LH=mu+6*sigma;
	for(i=0;i<n;i++)
	{
		double rs=gsl_vector_get(Res,i);
		if(rs>LH||rs<LL) gsl_vector_set(w,i,0.0);
	}	
	gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);
	gsl_multifit_linear_free(work);
	for(i=0;i<p;i++) par[i]=gsl_vector_get(c,i);
	double bReturn=gsl_RSquare(y,Res,w);
	gsl_vector_free(c);
	gsl_vector_free(w);
	gsl_vector_free(Res);
	gsl_matrix_free(cov);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	return bReturn;
}

int gsl_robustfit(gsl_matrix *X,gsl_vector *y, gsl_vector *c,double stats[3])
{
	gsl_matrix *cov;
	gsl_vector *w,*h,*Res;
	int p=X->size2;
	if(p<=1) return 0;
	int n=X->size1;
	double chisq;
	w=gsl_vector_alloc(n);
	h=gsl_vector_alloc(n);
	Res=gsl_vector_alloc(n);
	cov=gsl_matrix_alloc(p,p);
	int i;
	for(i=0;i<n;i++)
	{
		gsl_vector_set(w,i,1.0);
	}
	gsl_adfactor_get(X,h);
	gsl_multifit_linear_workspace * work= gsl_multifit_linear_alloc(n,p);	
	gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);	
    gsl_vector *c0;
	c0=gsl_vector_alloc(p);
	gsl_vector_set_zero(c0);	
	double max_cerr=1;
	double max_abs_c=0;
	double S_EPS=sqrt(EPS);
	double el,el1;
	int iter=0;
	while(max_cerr>S_EPS*max_abs_c)
	{
		gsl_res_get(X,y,c,Res);
		gsl_update_w(w,Res,h,8,p-1);
		gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);
		max_cerr=1;	
		max_abs_c=0;
		for(i=0;i<p;i++)
		{
			el=fabs(gsl_vector_get(c0,i)-gsl_vector_get(c,i));
			if(max_cerr>el)max_cerr=el;
			el=fabs(gsl_vector_get(c0,i));
			el1=fabs(gsl_vector_get(c,i));
			el=el>el1?el:el1;
			if(max_abs_c<el) max_abs_c=el;
		}
		gsl_vector_memcpy(c0,c);
		iter++;
		if(iter>MX_ITER_LIM) break;		
	}
	stats[0]=gsl_w_mean(Res,w);
	stats[1]=gsl_w_std(Res,w,stats[0]);
	stats[2]=gsl_RSquare(y,Res,w);
	gsl_multifit_linear_free(work);
	gsl_vector_free(c0);
	gsl_vector_free(h);
	gsl_vector_free(w);
	gsl_vector_free(Res);
	gsl_matrix_free(cov);
	return iter;
}

int gsl_robustfit(gsl_matrix *X,gsl_vector *y, gsl_vector *c)
{
	gsl_matrix *cov;
	gsl_vector *w,*h,*Res;
	int p=X->size2;
	if(p<=1) return 0;
	int n=X->size1;
	double chisq;
	w=gsl_vector_alloc(n);
	h=gsl_vector_alloc(n);
	Res=gsl_vector_alloc(n);
	cov=gsl_matrix_alloc(p,p);
	int i;
	for(i=0;i<n;i++)
	{
		gsl_vector_set(w,i,1.0);
	}
	gsl_adfactor_get(X,h);
	gsl_multifit_linear_workspace * work= gsl_multifit_linear_alloc(n,p);	
	gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);	
    gsl_vector *c0;
	c0=gsl_vector_alloc(p);
	gsl_vector_set_zero(c0);	
	double max_cerr=1;
	double max_abs_c=0;
	double S_EPS=sqrt(EPS);
	double el,el1;
	int count=0;
	int bReturn=n;
	while(max_cerr>S_EPS*max_abs_c)
	{
		gsl_res_get(X,y,c,Res);
		gsl_update_w(w,Res,h,8,p-1);
		gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);
		max_cerr=1;	
		max_abs_c=0;
		for(i=0;i<p;i++)
		{
			el=fabs(gsl_vector_get(c0,i)-gsl_vector_get(c,i));
			if(max_cerr>el)max_cerr=el;
			el=fabs(gsl_vector_get(c0,i));
			el1=fabs(gsl_vector_get(c,i));
			el=el>el1?el:el1;
			if(max_abs_c<el) max_abs_c=el;
		}
		gsl_vector_memcpy(c0,c);
		count++;
		if(count>100) 
		{
			bReturn=0;
			break;
		}
	}
	gsl_multifit_linear_free(work);
	gsl_vector_free(c0);
	gsl_vector_free(h);
	gsl_vector_free(w);
	gsl_vector_free(Res);
	gsl_matrix_free(cov);
	return bReturn;
}
////////////////////

double gsl_GetPA(double *x,double *y,int n)
{
	if(n<3) return 0;
	gsl_interp_accel *acc=gsl_interp_accel_alloc();
	gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline,n);
	gsl_spline_init(spline,x,y,n);
	double w=(x[n-1]-x[0])/500;
	double xi=x[0];
	double *y1;
	y1=new double[500];
	int i;
	for(i=0;i<500;i++)
	{		
		y1[i]=gsl_spline_eval(spline,xi,acc);
		xi+=w;
	}
	double PA=0;
	for(i=1;i<500;i++) PA+=(y1[i]+y1[i-1])*w/2;
	delete []y1;
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return PA;
}

double gsl_GetPCA(double *li,double *hi,int n)
{
	if(n<3) return 0;
	gsl_matrix * A;
	gsl_matrix * V;
	gsl_vector * S;
	gsl_vector * work;
	A=gsl_matrix_alloc(n,2);
	V=gsl_matrix_alloc(2,2);
	S=gsl_vector_alloc(2);
	work=gsl_vector_alloc(2);
	for(int i=0;i<n;i++)
	{
		gsl_matrix_set(A,i,0,li[i]);
		gsl_matrix_set(A,i,1,hi[i]);
	}
	gsl_linalg_SV_decomp(A,V,S,work);
	double P11=gsl_matrix_get(V,0,0);
	double P21=gsl_matrix_get(V,1,0);
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	if(fabs(P11)>DEPS) return P21/P11;
	else return 0;
}
//the size of x, y must great than n and nnew alwalys
bool gsl_PSmooth(double *x,double *y,int n,int nnew)
{
	if(n<=1) return false;
	if(nnew<n) return false;
	gsl_interp_accel *acc=gsl_interp_accel_alloc();
	gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline,n);
	gsl_spline_init(spline,x,y,n);
	double w=(x[n-1]-x[0])/nnew;
	double xi=x[0];	
	int i;
	for(i=0;i<nnew;i++)
	{		
		y[i]=gsl_spline_eval(spline,xi,acc);
		x[i]=xi;
		xi+=w;
	}
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return true;
}


double EM_Gauss_Fit(FIT_DATA *d,double beta[4],double L,double H)
{
	int Max_step=10;
	double Eps=1;
	double betaOLD[4];
	double r=Gauss_Fit(d,betaOLD);
	FIT_DATA dTmp;
	dTmp.x=new double[100];
	dTmp.y=new double[100];
	dTmp.n=100;
	double Step=(H-L)/100;
	int Iter=0;
	while(Eps>1e-4)
	{
		double xi=L;
		for(int i=0;i<100;i++)
		{
			xi+=Step;
			dTmp.x[i]=xi;
			dTmp.y[i]=gsl_gauss(dTmp.x[i],betaOLD);
		}
		r=Gauss_Fit(&dTmp,beta);
		Eps=0;
		for(int i=0;i<4;i++)
		{
			Eps+=fabs(betaOLD[i]-beta[i]);
			betaOLD[i]=beta[i];
		}
		Iter++;
		if(Iter>=Max_step) break;
	}
	return r;
}

//x, and y is the same size as default, not check here.
double Gauss_Fit(vector<double> &x,vector<double> &y,vector<double> &beta)
{
	vector<double> xbackup;
	vector<double> ybackup;
	int i,n=x.size();
	for(i=0;i<n;i++)
	{
		xbackup.push_back(x[i]);
		ybackup.push_back(y[i]);
	}	
	if(!Initial_cNew(xbackup,ybackup,beta)) return 0;
	OP op;
	op.DerivStep=exp(log(DEPS)/3);//这个参数增大会导致拟合失败
	op.Disp=2;
	op.MaxIter=500;//重新测试增加，有效
	op.TolFun=1e-8;
	op.TolX=1e-8;
	double *xTmp,*yTmp,*betaTmp;
	n=xbackup.size();
	if(n<=0) return 0;
	xTmp=new double[n];
	yTmp=new double[n];
	int nP=beta.size();
	betaTmp=new double[nP];
	for(i=0;i<n;i++) 
	{
		xTmp[i]=xbackup[i];
		yTmp[i]=ybackup[i];
	}
	for(i=0;i<nP;i++) betaTmp[i]=beta[i];
	gsl_nlfitNew(xTmp,yTmp,n,betaTmp,nP,&Gauss_F,&op);
	if(op.r>0.9)
	{
		for(i=0;i<nP;i++) beta[i]=betaTmp[i];
	}
	else if(nP==4&&n>4)//这个改进有效
	{	
		xbackup.clear();
		ybackup.clear();
		n=x.size();
		for(i=0;i<n;i++)
		{
			xbackup.push_back(x[i]);
			ybackup.push_back(y[i]);
		}

		Re_initial(xbackup,ybackup,beta);
		n=xbackup.size();
		if(n>0)
		{
			delete []betaTmp;
			delete []xTmp;
			delete []yTmp;

			betaTmp=new double[8];
			xTmp=new double[n];
			yTmp=new double[n];
			for(i=0;i<n;i++) 
			{
				xTmp[i]=xbackup[i];
				yTmp[i]=ybackup[i];
			}
			for(i=0;i<8;i++) betaTmp[i]=beta[i];
			gsl_nlfitNew(xTmp,yTmp,n,betaTmp,8,&Gauss_F,&op);
			if(op.r>0.9) for(i=0;i<8;i++) beta[i]=betaTmp[i];
		}
	}
	//如果两种拟合方案都失败，不更新beta，这样beta的值就是初始化的值
	delete []xTmp;
	delete []yTmp;
	delete []betaTmp;
	return op.r;
}

//find peaks in this cluster
	//          /|\      /|\       /|\      /|\   
	//         / | \    / | \     / | \    / | \                        
	//        /  |  \  /  |  \   /  |  \  /  |  \                       
	//       /   |   \/   |   \ /   |   \/   |   \                      
	//----------------------------------------------    
    //may be in the decrease acesent
	//          |\      /|\       /|\      /|   
	//          | \    / | \     / | \    / |                        
	//          |  \  /  |  \   /  |  \  /  |                       
	//          |   \/   |   \ /   |   \/   |                      
	//----------------------------------------------     

//x, and y is the same size as default, not check here.
//the structure of beta:
//int base_idx=4*i;
//    A=beta[base_idx];
//    mu=beta[base_idx+1];
//    sigma=beat[base_idx+2];
//    baseline=beta[base_idx+3];

#define  KPS  1.6652//sqrt(log(16.0))/2
bool Initial_c(vector<double> &x,vector<double> &y,vector<double> &beta)
{
	int i,n=x.size();
	if(n<=3)
	{
		double mu=0;
		double sigma=0;
		double A=0;
		for(i=0;i<n;i++)
		{
			mu+=x[i]*y[i];
			sigma+=y[i];
			if(A<y[i]) A=y[i];
		}	
		mu/=sigma;
		sigma=0.002;
		beta.push_back(A);
		beta.push_back(mu);
		beta.push_back(sigma);
		beta.push_back(0);
		return false;
	}	
	i=1;
	double MaxPK=0;
	vector<double> tmpBeta;
	while(1)
	{
		for(;i<n;i++)
		{
			if(y[i]<=y[i-1]) break;
		}	
		int pk=i-1;
		for(;i<n;i++)
		{
			if(y[i]>y[i-1]) break;
		}	
		int end=i-1;
		double sigma=(x[end]-x[pk])/KPS;
		if(sigma<0.0015) sigma=0.0015;//a default value according the experiment data
		//if(end<pk+2) sigma*=1.5;//add for little data points
		tmpBeta.push_back(y[pk]);
		tmpBeta.push_back(x[pk]);
		tmpBeta.push_back(sigma);
		tmpBeta.push_back(0.01*y[pk]); //base line estimation
		if(MaxPK<y[pk]) MaxPK=y[pk];
		if(i==n) break;
	}
	//truck some too low peaks
	int num=tmpBeta.size()/4;
	for(i=0;i<num;i++)
	{
		int baseidx=4*i;
		if(tmpBeta[baseidx]<0.05*MaxPK) continue;
		beta.push_back(tmpBeta[baseidx]);
		beta.push_back(tmpBeta[baseidx+1]);
		beta.push_back(tmpBeta[baseidx+2]);
		beta.push_back(tmpBeta[baseidx+3]);
	}
	
	//if(n<=3)
	//{
	//	beta[2]*=2;
	//}
	//please check the last peak, case special, but include in the workflow

	//check for not full peaks and interpration the data
    //method: a standard template was used to finish this process
	//
	vector<double> tmpX,tmpY;
	if(beta.size()<=0) return false;
	//extend the left
    double dm=x[1]-x[0];
	double xt=x[0]-dm;
	double BTtmp[4];
	for(i=0;i<4;i++) BTtmp[i]=beta[i];
	double yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	while(yt>0.01*(BTtmp[0]-BTtmp[3]))
	{
		tmpX.push_back(xt);
		tmpY.push_back(yt);
		xt-=dm;
		if(xt<BTtmp[1]-0.01) break;
		yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	}
	vector<double> swap;
	for(i=0;i<n;i++) swap.push_back(x[i]);
	x.clear();
	int nE=tmpX.size();
	for(i=nE-1;i>=0;i--)
	{
		x.push_back(tmpX[i]);
	}
	
	for(i=0;i<n;i++) x.push_back(swap[i]);

	swap.clear();
	for(i=0;i<n;i++) swap.push_back(y[i]);
	y.clear();
	nE=tmpY.size();
	for(i=nE-1;i>=0;i--)
	{
		y.push_back(tmpY[i]);
	}
	
	for(i=0;i<n;i++) y.push_back(swap[i]);

	x.push_back(xt);
	y.push_back(0);

	dm=x[n-1]-x[n-2];
	xt=x[n-1]+dm;
	int baseidx=beta.size()-4;
	for(i=0;i<4;i++) BTtmp[i]=beta[i+baseidx];
	yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	while(yt>0.01*(BTtmp[0]-BTtmp[3]))
	{
		x.push_back(xt);
		y.push_back(yt);
		xt+=dm;
		if(xt>BTtmp[1]+0.01) break;
		yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	}
	x.push_back(xt);
	y.push_back(0);
	return true;
}


bool Initial_cNew(vector<double> &x,vector<double> &y,vector<double> &beta)
{
	int i,n=x.size();
	if(n<=3)
	{
		double mu=0;
		double sigma=0;
		double A=0;
		for(i=0;i<n;i++)
		{
			mu+=x[i]*y[i];
			sigma+=y[i];
			if(A<y[i]) A=y[i];
		}	
		mu/=sigma;
		sigma=0.003;
		beta.push_back(A);
		beta.push_back(mu);
		beta.push_back(sigma);
		beta.push_back(0);
		return false;
	}	
	i=1;
	double MaxPK=0;
	vector<double> tmpBeta;
	while(1)
	{
		for(;i<n;i++)
		{
			if(y[i]<=y[i-1]) break;
		}	
		int pk=i-1;
		for(;i<n;i++)
		{
			if(y[i]>y[i-1]) break;
		}	
		int end=i-1;
		double sigma=(x[end]-x[pk])/KPS;
		if(sigma<0.002) sigma=0.002;//a default value according the experiment data
		//if(end<pk+2) sigma*=1.5;//add for little data points
		tmpBeta.push_back(y[pk]);
		tmpBeta.push_back(x[pk]);
		tmpBeta.push_back(sigma);
		tmpBeta.push_back(0.01*y[pk]); //base line estimation
		if(MaxPK<y[pk]) MaxPK=y[pk];
		if(i==n) break;
	}
	//truck some too low peaks
	int num=tmpBeta.size()/4;
	for(i=0;i<num;i++)
	{
		int baseidx=4*i;
		if(tmpBeta[baseidx]<0.05*MaxPK) continue;
		beta.push_back(tmpBeta[baseidx]);
		beta.push_back(tmpBeta[baseidx+1]);
		beta.push_back(tmpBeta[baseidx+2]);
		beta.push_back(tmpBeta[baseidx+3]);
	}
	
	//if(n<=3)
	//{
	//	beta[2]*=2;
	//}
	//please check the last peak, case special, but include in the workflow

	//check for not full peaks and interpration the data
    //method: a standard template was used to finish this process
	//
	vector<double> tmpX,tmpY;
	if(beta.size()<=0) return false;
	//extend the left
    double dm=x[1]-x[0];
	if(dm<0.0005) dm=0.0005;
	double xt=x[0]-dm;
	double BTtmp[4];
	for(i=0;i<4;i++) BTtmp[i]=beta[i];
	double yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	while(yt>0.01*(BTtmp[0]-BTtmp[3]))
	{
		tmpX.push_back(xt);
		tmpY.push_back(yt);
		xt-=dm;
		if(xt<BTtmp[1]-0.01) break;
		yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	}
	vector<double> swap;
	for(i=0;i<n;i++) swap.push_back(x[i]);
	x.clear();
	int nE=tmpX.size();
	for(i=nE-1;i>=0;i--)
	{
		x.push_back(tmpX[i]);
	}
	
	for(i=0;i<n;i++) x.push_back(swap[i]);

	swap.clear();
	for(i=0;i<n;i++) swap.push_back(y[i]);
	y.clear();
	nE=tmpY.size();
	for(i=nE-1;i>=0;i--)
	{
		y.push_back(tmpY[i]);
	}
	
	for(i=0;i<n;i++) y.push_back(swap[i]);

	x.push_back(xt);
	y.push_back(0);

	dm=x[n-1]-x[n-2];
	if(dm<0.0005) dm=0.0005;
	xt=x[n-1]+dm;
	int baseidx=beta.size()-4;
	for(i=0;i<4;i++) BTtmp[i]=beta[i+baseidx];
	yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	while(yt>0.01*(BTtmp[0]-BTtmp[3]))
	{
		x.push_back(xt);
		y.push_back(yt);
		xt+=dm;
		if(xt>BTtmp[1]+0.01) break;
		yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	}
	x.push_back(xt);
	y.push_back(0);
	return true;
}

//the size of beta must be 4, x and y must be original data, not extend data
void Re_initial(vector<double> &x,vector<double> &y,vector<double> &beta)
{
	int i,n=x.size();
	i=0;
	double MaxDiff=0;
	for(int k=0;k<n;k++)
	{
		double yt=Gauss_F(x[k],beta);
		yt=y[k]-yt;
		if(yt>MaxDiff)
		{
			i=k;
			MaxDiff=yt;
		}
	}
	beta[2]/=1.5;
	beta.push_back(y[i]);
	beta.push_back(x[i]);
	beta.push_back(beta[2]/1.5);
	beta.push_back(0.01*y[i]);
	//sort the beta according the mu
	if(beta[1]>beta[5])
	{
		for(i=0;i<4;i++)
		{
			double swap=beta[i];
			beta[i]=beta[4+i];
			beta[4+i]=swap;
		}
	}
	//re extend the data again
	vector<double> tmpX,tmpY;
	//extend the left
    double dm=x[1]-x[0];
	if(dm<0.0005) dm=0.0005;
	double xt=x[0]-dm;
	double BTtmp[4];
	for(i=0;i<4;i++) BTtmp[i]=beta[i];
	double yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	while(yt>0.01*(BTtmp[0]-BTtmp[3]))
	{
		tmpX.push_back(xt);
		tmpY.push_back(yt);
		xt-=dm;
		if(xt<BTtmp[1]-0.01) break;
		yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	}
	vector<double> swap;
	for(i=0;i<n;i++) swap.push_back(x[i]);
	x.clear();
	int nE=tmpX.size();
	for(i=nE-1;i>=0;i--)
	{
		x.push_back(tmpX[i]);
	}
	
	for(i=0;i<n;i++) x.push_back(swap[i]);

	swap.clear();
	for(i=0;i<n;i++) swap.push_back(y[i]);
	y.clear();
	nE=tmpY.size();
	for(i=nE-1;i>=0;i--)
	{
		y.push_back(tmpY[i]);
	}
	
	for(i=0;i<n;i++) y.push_back(swap[i]);

	x.push_back(xt);
	y.push_back(0);

	dm=x[n-1]-x[n-2];
	if(dm<0.0005) dm=0.0005;
	xt=x[n-1]+dm;
	int baseidx=beta.size()-4;
	for(i=0;i<4;i++) BTtmp[i]=beta[i+baseidx];
	yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	while(yt>0.01*(BTtmp[0]-BTtmp[3]))
	{
		x.push_back(xt);
		y.push_back(yt);
		xt+=dm;
		if(xt>BTtmp[1]+0.01) break;
		yt=Gauss_F(xt,BTtmp)-BTtmp[3];
	}
	x.push_back(xt);
	y.push_back(0);	
}

//bool Initial_c(vector<double> &x,vector<double> &y,vector<double> &beta)
//{
//	int i,n=x.size();
//	if(n<=2) return false;
//	i=1;
//	int begin=0;
//	while(1)
//	{
//		for(;i<n;i++)
//		{
//			if(y[i]<=y[i-1]) break;
//		}	
//		int pk=i-1;
//		for(;i<n;i++)
//		{
//			if(y[i]>y[i-1]) break;
//		}	
//		//i--;
//		if(i==n) i--;
//		double sigma=(x[i]-x[begin])/KPS;
//		if(sigma<=0.0015) sigma=0.0015;//a default value according the experiment data	
//		//double sigma=0.015;
//		//sigma=0.024007943396438865;
//		beta.push_back(y[pk]);
//		beta.push_back(x[pk]);
//		beta.push_back(sigma);
//		beta.push_back(0.01*y[pk]); //base line estimation
//		if(i==n-1) break;
//		begin=i;
//	}
//	//please check the last peak, case special, but include in the workflow
//
//	//check for not full peaks and interpration the data
//    //method: a standard template was used to finish this process
//	//
//	/*vector<double> tmpX,tmpY;*/
//	int num=beta.size();
//	if(num<4) return false;
//	//extend the left
//    double dm=beta[2]*sqrt(log(100*beta[0]));
//	double xt=beta[1]-dm;	
//	//xt=473.19000244140625;debug
//	vector<double> swap;
//	for(i=0;i<n;i++) swap.push_back(x[i]);
//	x.clear();
//	x.push_back(xt);
//	for(i=0;i<n;i++) x.push_back(swap[i]);
//
//	swap.clear();
//	for(i=0;i<n;i++) swap.push_back(y[i]);
//	y.clear();
//	y.push_back(0);
//
//	for(i=0;i<n;i++) y.push_back(swap[i]);
//	
//	dm=beta[num-2]*sqrt(log(100*beta[num-4]));
//	xt=beta[num-3]+dm;
//	//xt=473.22998046875000;debug
//	x.push_back(xt);
//	y.push_back(0);	
//	return true;
//}
double Gauss_F(double x,double *beta,int num)
{
	double sum=0;
	for(int i=0;i<num;i++)
	{
		int base_idx=4*i;
		double A=beta[base_idx];
		double mu=beta[base_idx+1];
		double sigma=beta[base_idx+2];
		double baseline=beta[base_idx+3];
		double t=(x-mu)/sigma;
		sum+=A*exp(-t*t)+baseline;		
	}
	return sum;
}

double Gauss_F(double x, vector<double> &beta)
{
	double sum=0;
	int num=beta.size()/4;
	for(int i=0;i<num;i++)
	{
		int base_idx=4*i;
		double A=beta[base_idx];
		double mu=beta[base_idx+1];
		double sigma=beta[base_idx+2];
		double baseline=beta[base_idx+3];
		double t=(x-mu)/sigma;
		sum+=A*exp(-t*t)+baseline;		
	}
	return sum;
}
////////////////////////////////new fitting functions
double Gauss_F_new(double x,double *beta,int num)
{
	double sum=0;
	int idx=0;
	for(int i=0;i<num;i++)
	{	
		double A=beta[idx];
		idx++;
		double mu=beta[idx];
		idx++;
		double sigma=beta[idx];
		idx++;		
		double t=(x-mu)/sigma;
		sum+=A*exp(-t*t);		
	}
	sum+=beta[idx];
	return sum;
}

double Gauss_F_new(double x, vector<double> &beta)
{
	double sum=0;
	int num=beta.size()/3;
	int idx=0;
	for(int i=0;i<num;i++)
	{		
		double A=beta[idx];
		idx++;
		double mu=beta[idx];
		idx++;
		double sigma=beta[idx];	
		idx++;
		double t=(x-mu)/sigma;
		sum+=A*exp(-t*t);		
	}
	sum+=beta[idx];
	return sum;
}

int guass_f(const gsl_vector * par, void *data,gsl_vector * f)
{
	size_t n = ((struct guass_data *)data)->n;
	double *y = ((struct guass_data *)data)->y;
	double *weight = ((struct guass_data *) data)->weight;
	double *x=((struct guass_data *) data)->x;
	size_t fn=((struct guass_data *) data)->fn;
	//do not check fn, at least 1
	double *A,*mu,*sigma,base;
	A=new double[fn];
	mu=new double[fn];
	sigma=new double[fn];
	size_t i,k;
	k=0;
	for(i=0;i<fn;i++)
	{
		A[i]=gsl_vector_get (par,k);
		k++;
		mu[i]=gsl_vector_get (par,k);
		k++;
		sigma[i]=gsl_vector_get (par,k);
		k++;
	}
	double t;
	base=gsl_vector_get (par,k);
	for (i = 0;i<n;i++)
	{
		/* Model Yi = sum(K,A * exp-((Xi-muK)/sigmaK)^2 )+ base */
		double Yi=base;
		for(k=0;k<fn;k++)
		{
			t=(x[i]-mu[k])/sigma[k];
			Yi+=A[k]*exp(-t*t);
		}		
		gsl_vector_set (f,i,(Yi-y[i])/weight[i]);
	}
	delete []A;
	delete []mu;
	delete []sigma;
	return GSL_SUCCESS;
}

int guass_df(const gsl_vector * par, void *data,gsl_matrix * J)
{
	size_t n = ((struct guass_data *)data)->n;
	double *y = ((struct guass_data *)data)->y;
	double *weight = ((struct guass_data *) data)->weight;
	double *x=((struct guass_data *) data)->x;
	size_t fn=((struct guass_data *) data)->fn;
	//do not check fn, at least 1
	double *A,*mu,*sigma;
	A=new double[fn];
	mu=new double[fn];
	sigma=new double[fn];
	size_t i,k;
	k=0;
	for(i=0;i<fn;i++)
	{
		A[i]=gsl_vector_get (par,k);
		k++;
		mu[i]=gsl_vector_get (par,k);
		k++;
		sigma[i]=gsl_vector_get (par,k);
		k++;
	}
	double t;	
	for(i=0;i<n;i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i], */
		/* Yi = sum(K,A * exp-((Xi-muK)/sigmaK)^2 )+ base */
		/* and the xj are the parameters (A,mu,sigma,base) */
		double s = weight[i];
		for(k=0;k<fn;k++)
		{
			t=(x[i]-mu[k])/sigma[k];			
			double e = exp(-t*t);
			gsl_matrix_set (J, i, 0, e/s);
			e*=2*A[k]*t/sigma[k];
			gsl_matrix_set (J, i, 1,e/s);
			e*=t;
			gsl_matrix_set (J, i, 2,e/s);
		}
		gsl_matrix_set (J, i, fn*3,1/s);
	}
	delete []A;
	delete []mu;
	delete []sigma;
	return GSL_SUCCESS;
}

int guass_fdf (const gsl_vector * x, void *data,gsl_vector * f, gsl_matrix * J)
{
	guass_f (x,data,f);
	guass_df (x,data,J);
	return GSL_SUCCESS;
}


int gsl_nlfit_og(double *x,double *y,const size_t n,double *beta,size_t fn,gsl_nlfit_options *op)
{
	if(fn<=0||n<=0) return 0;
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	size_t i, iter = 0;	
	double *weight;
	weight=new double[n];	
	for(i=0;i<n;i++) weight[i]=1.0;
	struct guass_data d = {fn,n,x,y,weight};
	gsl_multifit_function_fdf f;	
	const size_t p = 3*fn+1;

	gsl_vector *par;
	par=gsl_vector_alloc(p);
	for(i=0;i<p;i++) gsl_vector_set(par,i,beta[i]);
	//gsl_vector_view par = gsl_vector_view_array(beta,p);
	f.f = &guass_f;
	f.df = &guass_df;
	f.fdf = &guass_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;	
	
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, n, p);
	gsl_multifit_fdfsolver_set(s, &f, par);
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);		
		if(status==GSL_SUCCESS)
		{
			if(iter>3)	break;
			else status=GSL_CONTINUE;
		}
		status = gsl_multifit_test_delta (s->dx,s->x,op->TolFun,op->TolX);
	}while(status == GSL_CONTINUE&&iter<op->MaxIter);
		
	if(status == GSL_SUCCESS)//not good keep the initial
	{		
		for(i=0;i<p;i++)beta[i]=gsl_vector_get(s->x,i);
	}	
	//gsl_vector *fc;
	//fc=gsl_vector_alloc(n);
	//guass_f(s->x,&d,fc);
	double sse=gsl_blas_dnrm2(s->f);
	sse=sse*sse;
	double sst=gsl_stats_sd(y,1,n);
	sst=sst*sst;
	op->r=1-sse/sst;	
	gsl_multifit_fdfsolver_free(s);
	gsl_vector_free(par);
	//gsl_vector_free(fc);
	delete []weight;
	return iter;
}


bool Initial_par(vector<double> &x,vector<double> &y,vector<double> &beta)
{
	int i,n=x.size();
	if(n<=3)
	{
		double mu=0;
		double sigma=0;
		double A=0;
		for(i=0;i<n;i++)
		{
			mu+=x[i]*y[i];
			sigma+=y[i];
			if(A<y[i]) A=y[i];
		}	
		mu/=sigma;
		sigma=0.0015;
		beta.push_back(A);
		beta.push_back(mu);
		beta.push_back(sigma);
		beta.push_back(0.001*A);
		return false;
	}	
	i=1;
	double MaxPK=0;
	vector<double> tmpBeta;
	while(1)
	{
		for(;i<n;i++)
		{
			if(y[i]<=y[i-1]) break;
		}	
		int pk=i-1;
		for(;i<n;i++)
		{
			if(y[i]>y[i-1]) break;
		}	
		int end=i-1;
		double sigma=(x[end]-x[pk])/KPS;
		if(sigma<0.0015) sigma=0.0015;//a default value according the experiment data
		//if(end<pk+2) sigma*=1.5;//add for little data points
		tmpBeta.push_back(y[pk]);
		tmpBeta.push_back(x[pk]);
		tmpBeta.push_back(sigma);	
		if(MaxPK<y[pk]) MaxPK=y[pk];
		if(i==n) break;
	}
	//truck some too low peaks
	int num=tmpBeta.size()/3;
	int idx=0;
	for(i=0;i<num;i++)
	{	
		if(tmpBeta[idx]<0.01*MaxPK) 
		{
			idx+=3;
			continue;
		}
		//idx++;
		beta.push_back(tmpBeta[idx]);
		idx++;
		beta.push_back(tmpBeta[idx]);
		idx++;
		beta.push_back(tmpBeta[idx]);	
		idx++;
	}

	beta.push_back(MaxPK*0.001);//the overlall baseline	

	//please check the last peak, case special, but include in the workflow

	//check for not full peaks and interpration the data
    //method: a standard template was used to finish this process
	//
	vector<double> tmpX,tmpY;
	if(beta.size()<=0) return false;
	//extend the left
    double dm=x[1]-x[0];
	double xt=x[0]-dm;
	double BTtmp[4];
	for(i=0;i<3;i++) BTtmp[i]=beta[i];
	BTtmp[3]=0;//baseline estimate
	double yt=Gauss_F_new(xt,BTtmp,1);
	while(yt>0.01*BTtmp[0])
	{
		tmpX.push_back(xt);
		tmpY.push_back(yt);
		xt-=dm;
		if(xt<BTtmp[1]-0.01) break;
		yt=Gauss_F_new(xt,BTtmp,1);
	}
	vector<double> swap;
	for(i=0;i<n;i++) swap.push_back(x[i]);
	x.clear();
	int nE=tmpX.size();
	for(i=nE-1;i>=0;i--)
	{
		x.push_back(tmpX[i]);
	}
	
	for(i=0;i<n;i++) x.push_back(swap[i]);

	swap.clear();
	for(i=0;i<n;i++) swap.push_back(y[i]);
	y.clear();
	nE=tmpY.size();
	for(i=nE-1;i>=0;i--)
	{
		y.push_back(tmpY[i]);
	}
	
	for(i=0;i<n;i++) y.push_back(swap[i]);

	x.push_back(xt);
	y.push_back(0);

	dm=x[n-1]-x[n-2];
	xt=x[n-1]+dm;
	int baseidx=beta.size()-4;
	for(i=0;i<4;i++) BTtmp[i]=beta[i+baseidx];
	yt=Gauss_F_new(xt,BTtmp,1);
	while(yt>0.01*BTtmp[0])
	{
		x.push_back(xt);
		y.push_back(yt);
		xt+=dm;
		if(xt>BTtmp[1]+0.01) break;
		yt=Gauss_F_new(xt,BTtmp,1);
	}
	x.push_back(xt);
	y.push_back(0);
	return true;
}

//the size of beta must be 4, x and y must be original data, not extend data
void Re_initial_par(vector<double> &x,vector<double> &y,vector<double> &beta)
{
	int i,n=x.size();
	i=0;
	double MaxDiff=0;
	for(int k=0;k<n;k++)
	{
		double yt=Gauss_F_new(x[k],beta);
		yt=y[k]-yt;
		if(yt>MaxDiff)
		{
			i=k;
			MaxDiff=yt;
		}
	}
	beta[2]/=1.5;
	beta.push_back(y[i]);
	beta.push_back(x[i]);
	beta.push_back(beta[2]);	
	//sort the beta according the mu
	if(beta[1]>beta[4])
	{
		for(i=0;i<3;i++)
		{
			double swap=beta[i];
			beta[i]=beta[3+i];
			beta[3+i]=swap;
		}
	}
	//re extend the data again
	vector<double> tmpX,tmpY;
	//extend the left
    double dm=x[1]-x[0];
	double xt=x[0]-dm;
	double BTtmp[4];
	for(i=0;i<3;i++) BTtmp[i]=beta[i];
	BTtmp[3]=0;
	double yt=Gauss_F_new(xt,BTtmp,1);
	double LowLim=BTtmp[0];
	if(LowLim<0.01) LowLim=0.01;
	while(yt>LowLim)
	{
		tmpX.push_back(xt);
		tmpY.push_back(yt);
		xt-=dm;
		if(xt<BTtmp[1]-0.01) break;
		yt=Gauss_F_new(xt,BTtmp,1);
	}
	vector<double> swap;
	for(i=0;i<n;i++) swap.push_back(x[i]);
	x.clear();
	int nE=tmpX.size();
	for(i=nE-1;i>=0;i--)
	{
		x.push_back(tmpX[i]);
	}
	
	for(i=0;i<n;i++) x.push_back(swap[i]);

	swap.clear();
	for(i=0;i<n;i++) swap.push_back(y[i]);
	y.clear();
	nE=tmpY.size();
	for(i=nE-1;i>=0;i--)
	{
		y.push_back(tmpY[i]);
	}
	
	for(i=0;i<n;i++) y.push_back(swap[i]);

	x.push_back(xt);
	y.push_back(0);

	dm=x[n-1]-x[n-2];
	xt=x[n-1]+dm;
	int baseidx=beta.size()-4;
	for(i=0;i<3;i++) BTtmp[i]=beta[i+baseidx];
	yt=Gauss_F_new(xt,BTtmp,1);
	while(yt>0.01*BTtmp[0])
	{
		x.push_back(xt);
		y.push_back(yt);
		xt+=dm;
		if(xt>BTtmp[1]+0.01) break;
		yt=Gauss_F_new(xt,BTtmp,1);
	}
	x.push_back(xt);
	y.push_back(0);	
}


//x, and y is the same size as default, not check here.
double Gauss_Fit_new(vector<double> &x,vector<double> &y,vector<double> &beta)
{
	vector<double> xbackup;
	vector<double> ybackup;
	int i,n=x.size();
	for(i=0;i<n;i++)
	{
		xbackup.push_back(x[i]);
		ybackup.push_back(y[i]);
	}	
	if(!Initial_par(x,y,beta)) return 0;
	OP op;
	op.MaxIter=200;//重新测试增加，有效,500
	op.TolFun=1e-6;
	op.TolX=1e-6;
	double *xTmp,*yTmp,*betaTmp;
	n=x.size();
	if(n<=0) return 0;
	xTmp=new double[n];
	yTmp=new double[n];
	int nP=beta.size();
	betaTmp=new double[nP];
	for(i=0;i<n;i++) 
	{
		xTmp[i]=x[i];
		yTmp[i]=y[i];
	}
	for(i=0;i<nP;i++) betaTmp[i]=beta[i];
	gsl_nlfit_og(xTmp,yTmp,n,betaTmp,nP/3,&op);	
	if(op.r<0.9&&nP==4&&n>4)//这个改进有效
	{	
		Re_initial_par(xbackup,ybackup,beta);
		n=xbackup.size();
		if(n>0)
		{
			delete []betaTmp;
			delete []xTmp;
			delete []yTmp;

			betaTmp=new double[8];
			xTmp=new double[n];
			yTmp=new double[n];
			for(i=0;i<n;i++) 
			{
				xTmp[i]=xbackup[i];
				yTmp[i]=ybackup[i];
			}
			for(i=0;i<7;i++) betaTmp[i]=beta[i];
			gsl_nlfit_og(xTmp,yTmp,n,betaTmp,2,&op);
			if(op.r>0.9) for(i=0;i<7;i++) beta[i]=betaTmp[i];//modifed on 2010.11.2			
		}
	}
	else 
	{
		for(i=0;i<nP;i++) beta[i]=betaTmp[i];//move to here to aviod the abnormal value of betaTmp
	}
	delete []xTmp;
	delete []yTmp;
	delete []betaTmp;
	return op.r;
}

int weight_centriod(vector<double> &x,vector<double> &y)
{		
	int i,n=x.size();
	if(n<=2) return n;
	size_t begin=0;
	size_t end=1;
	double tmz,tint;
	vector<double> mz_c;
	vector<double> int_c;
	while(begin<n)
	{
		for(i=begin+1;i<n;i++)
		{
			if(y[i]<y[i-1]) break;
			if(x[i]>x[i-1]+0.05) break;
		}
		for(;i<n;i++)
		{
			if(y[i]>y[i-1]) break;
			if(x[i]>x[i-1]+0.05) break;
		}
		end=i;
		//if(end==n) end--;
		if((end<n-1)&(end>1+begin))		
		{
			double d1=y[end]-y[end-1];
			double d2=y[end-2]-y[end-1];
			if(d2>d1) end--;
		}
		tmz=0;
		tint=0;
		for(i=begin;i<end;i++)
		{
			tmz+=x[i]*y[i];
			tint+=y[i];
		}
		tmz/=tint;
		double dmax=x[begin]-tmz;
		double tf=x[end-1]-tmz;
		dmax=dmax>tf?dmax:tf;
		for(i=begin;i<end;i++)
		{
			tf=(x[i]-tmz)/dmax;
			tint+=y[i]*exp(-tf*tf);
		}
		mz_c.push_back(tmz);
		int_c.push_back(tint);		
		begin=end;
	}
	x.clear();
	y.clear();
	n=mz_c.size();
	for(i=0;i<n;i++)
	{
		x.push_back(mz_c[i]);
		y.push_back(int_c[i]);
	}	
	return n;
}

double gsl_stats_correlation (const double data1[], const size_t stride1,const double data2[], const size_t stride2,const size_t n)
{
  size_t i;
  long double sum_xsq = 0.0;
  long double sum_ysq = 0.0;
  long double sum_cross = 0.0;
  long double ratio;
  long double delta_x, delta_y;
  long double mean_x, mean_y;
  long double r;

  /*
   * Compute:
   * sum_xsq = Sum [ (x_i - mu_x)^2 ],
   * sum_ysq = Sum [ (y_i - mu_y)^2 ] and
   * sum_cross = Sum [ (x_i - mu_x) * (y_i - mu_y) ]
   * using the above relation from Welford's paper
   */

  mean_x = data1[0 * stride1];
  mean_y = data2[0 * stride2];

  for (i = 1; i < n; ++i)
    {
      ratio = i / (i + 1.0);
      delta_x = data1[i * stride1] - mean_x;
      delta_y = data2[i * stride2] - mean_y;
      sum_xsq += delta_x * delta_x * ratio;
      sum_ysq += delta_y * delta_y * ratio;
      sum_cross += delta_x * delta_y * ratio;
      mean_x += delta_x / (i + 1.0);
      mean_y += delta_y / (i + 1.0);
    }

  r = sum_cross / (sqrt(sum_xsq) * sqrt(sum_ysq));

  return r;
}


double mygsl_stats_mean(vector<double> &ppme)
{
	double *x;
	size_t i,n=ppme.size();
	if(n<=0) return 0;
	x=new double[n];
	for(i=0;i<n;i++) x[i]=ppme[i];
	double m=gsl_stats_mean(x,1,n);
	delete []x;
	return m;

}

double mygsl_stats_std_m(vector<double> &ppme,double m)
{
	double *x;
	size_t i,n=ppme.size();
	if(n<=0) return 0;
	x=new double[n];
	for(i=0;i<n;i++) x[i]=ppme[i];
	double sd=gsl_stats_sd_m(x,1,n,m);
	delete []x;
	return sd;
}

void mygsl_stats_stdm(vector<double> &ppme,double m[2])
{
	double *x;
	size_t i,n=ppme.size();
	m[1]=0;
	m[2]=0;
	if(n<=0) return;
	x=new double[n];
	for(i=0;i<n;i++) x[i]=ppme[i];
	m[1]=gsl_stats_mean(x,1,n);
	m[2]=gsl_stats_sd_m(x,1,n,m[1]);
	delete []x;	
}

void mygsl_stats_stdm(vector<double> &ppme,double *m,double *sd)
{
	double *x;
	size_t i,n=ppme.size();
	m[1]=0;
	m[2]=0;
	if(n<=0) return;
	x=new double[n];
	for(i=0;i<n;i++) x[i]=ppme[i];
	*m=gsl_stats_mean(x,1,n);
	*sd=gsl_stats_sd_m(x,1,n,m[1]);
	delete []x;	
}

