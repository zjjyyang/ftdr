#include "StdAfx.h"
#include "IsoGroupD.h"

#include "gsl/gsl_multimin.h"
#include "GaussFit.h"
#include "OTrace.h"
//	
sPeak::sPeak()
{
	mz=0;
	dInt=0;
	base=0;
	noise=0;
	IsUsed=false;
	Global_IDX=-1;	
	count=0;
}

sPeak::~sPeak()
{

}

sPeak::sPeak(const sPeak &r)
{
	mz=r.mz;
	dInt=r.dInt;	
	base=r.base;
	noise=r.noise;
	IsUsed=r.IsUsed;
	Global_IDX=r.Global_IDX;	
	size_t sg=r.LocalList.size();
	for(size_t i=0;i<sg;i++) LocalList.push_back(r.LocalList[i]);
	count=r.count;
}

bool sPeak::IsInLocalList(int idx)
{
	size_t sg=LocalList.size();
	for(size_t i=0;i<sg;i++) 
	{
		if(idx==LocalList[i].UnitIdx) return true;
	}
	return false;
}

sPeak &sPeak::operator=(const sPeak &r)
{
	mz=r.mz;
	dInt=r.dInt;
	base=r.base;
	noise=r.noise;
	IsUsed=r.IsUsed;
	Global_IDX=r.Global_IDX;
	LocalList.clear();
	size_t sg=r.LocalList.size();
	for(size_t i=0;i<sg;i++) LocalList.push_back(r.LocalList[i]);
	count=r.count;
	return *this;
}

sPeak &sPeak::operator=(mPeak &r)
{
	mz=r.dMass;
	dInt=r.dIntensity;
	base=r.baseLine;
	noise=r.noise;
	IsUsed=false;
	Global_IDX=-1;
	LocalList.clear();	
	count=0;
	return *this;
}

sPeak::sPeak(const mPeak &r)
{
	mz=r.dMass;
	dInt=r.dIntensity;
	base=r.baseLine;
	noise=r.noise;
	IsUsed=false;
	Global_IDX=-1;
	LocalList.clear();	
	count=0;
}

refPeak::refPeak()
{
	mz=0;
	dInt=0;
	clusterid=-1;
	localid=-1;
	IsoNum=0;
	S_N=0;
	rank=-1;
	weight=1.0;
}

refPeak::~refPeak()
{

}

refPeak::refPeak(const refPeak &cp)
{
	mz=cp.mz;
	dInt=cp.dInt;
	clusterid=cp.clusterid;
	localid=cp.localid;
	IsoNum=cp.IsoNum;
	S_N=cp.S_N;
	rank=cp.rank;
	weight=cp.weight;
	
}

refPeak &refPeak::operator=(const refPeak &cp)
{
	mz=cp.mz;
	dInt=cp.dInt;	
	clusterid=cp.clusterid;
	localid=cp.localid;
	IsoNum=cp.IsoNum;
	S_N=cp.S_N;
	rank=cp.rank;
	weight=cp.weight;
	return *this;
}

refPeak &refPeak::operator=(const sPeak &cp)
{
	mz=cp.mz;
	dInt=cp.dInt;	
	S_N=cp.base+cp.noise;	
	clusterid=-1;
	localid=-1;
	if(S_N>1e-6) S_N=dInt/S_N;
	else S_N=0;
	rank=-1;
	weight=1.0;
	return *this;
}


FitVector::FitVector()
{

}

FitVector::~FitVector()
{

}

FitVector::FitVector(const FitVector &r)
{
	Exp_Int=r.Exp_Int;
	size_t i,gs=r.IDX.size();
	//IDX.clear();
	for(i=0;i<gs;i++) IDX.push_back(r.IDX.at(i));
	for(i=0;i<gs;i++) IsoTL.push_back(r.IsoTL.at(i));	
}

FitVector &FitVector::operator=(const FitVector &r)
{
	Exp_Int=r.Exp_Int;
	size_t i,gs=r.IDX.size();
	IDX.clear();
	IsoTL.clear();
	for(i=0;i<gs;i++) IDX.push_back(r.IDX.at(i));
	for(i=0;i<gs;i++) IsoTL.push_back(r.IsoTL.at(i));	
	return *this;
}

int FitVector:: FindIdx(int idx)
{
	size_t i,gs=IDX.size();
	for(i=0;i<gs;i++)
		if(IDX[i]==idx) return i;
	return -1;
}

//	
IsoUnit::IsoUnit()
{
	charge=CH_UNKNOWN;
	weight=1.0;
}

IsoUnit::~IsoUnit()
{

}

IsoUnit::IsoUnit(const IsoUnit &cp)
{
	charge=cp.charge;
	weight=cp.weight;
	size_t i,sg=cp.peak_idx.size();
	for(i=0;i<sg;i++) peak_idx.push_back(cp.peak_idx[i]);
}

IsoUnit &IsoUnit::operator=(const IsoUnit &cp)
{
	peak_idx.clear();
	charge=cp.charge;
	weight=cp.weight;
	size_t i,sg=cp.peak_idx.size();
	for(i=0;i<sg;i++) peak_idx.push_back(cp.peak_idx[i]);
	return *this;
}

bool IsoUnit::CheckFactor(IsoUnit &cs)
{
	if(charge==cs.charge) return false;
	if(cs.charge%charge!=0) return false;	
	size_t i,sg=cs.peak_idx.size();
	size_t k,sg1=peak_idx.size();
	if(sg<=sg1) return false;
	for(i=0;i<sg1;i++)
	{
		for(k=0;k<sg;k++)
		{
			if(peak_idx[i]==cs.peak_idx[k]) break;
		}
		if(k==sg) return false;
	}
	return true;
}

bool IsoUnit::CheckCover(IsoUnit &cs)
{
	if(charge!=cs.charge) return false;	
	size_t i,sg=cs.peak_idx.size();
	size_t k,sg1=peak_idx.size();
	if(sg>=sg1) return false;
	for(i=0;i<sg;i++)
	{
		for(k=0;k<sg1;k++)
		{
			if(peak_idx[i]==cs.peak_idx[k]) break;
		}
		if(k==sg1) return false;
	}
	return true;
}

void IsoUnit::AdjustCH(double mass)
{
	size_t i,sg=peak_idx.size();
	if(sg<=0) return;
	for(i=0;i<sg;i++)
	{
		if(GetIsoDis(mass,i+1)<1e-4) break;
	}	
	for(;sg>i;sg--) peak_idx.pop_back();//find a bug on 2011.11.29,remove the =
}

IsoLocal::IsoLocal()
{
	UnitIdx=-1;
	UnitLocal=-1;
}

IsoLocal::~IsoLocal()
{

}

IsoLocal::IsoLocal(const IsoLocal &cp)
{
	UnitIdx=cp.UnitIdx;
	UnitLocal=cp.UnitLocal;

}

IsoLocal& IsoLocal::operator =(const IsoLocal &cp)
{
	UnitIdx=cp.UnitIdx;
	UnitLocal=cp.UnitLocal;
	return *this;

}

bool IsoGroup::InitialPars(vector<refPeak> &ACIs,fParList &PL)
{
	PL.clear();
	FitVector ft;
	size_t k,i,sNum=UList.size();
	if(sNum<=0) return false;
	size_t pnum=gPeaks.size();
	if(pnum<=0) return false;
	for(i=0;i<pnum;i++)
	{
		ft.Exp_Int=gPeaks[i].dInt;
		PL.push_back(ft);
	}

	ACIs.clear();
	refPeak rt;

	int tid;
	int curAIdx=0;
	size_t IS;
	double IsoDisT[MAX_ISO];
	//IG.SortUList();


	//add the charge determined overlap at first
	for(i=0;i<sNum;i++)
	{
		IS=UList[i].peak_idx.size();
		if(IS<=0) continue;
		if(UList[i].weight<0.9) continue;
		tid=UList[i].peak_idx[0];
		PL[tid].IDX.push_back(curAIdx);
		rt=gPeaks[tid];
		ACIs.push_back(rt);
		size_t ts=ACIs.size()-1;
		ACIs[ts].clusterid=i;
		ACIs[ts].localid=0;	
		ACIs[ts].IsoNum=IS;
		ACIs[ts].weight=UList[i].weight;
		double mass=gPeaks[tid].mz*UList[i].charge;//process the unknown charge
	
		GetIsoDis(mass,IsoDisT);//if charge is 0, this function will only return IsoDisT[0]=1.0;, IS=1;
		PL[tid].IsoTL.push_back(IsoDisT[0]);
		for(k=1;k<IS;k++)
		{
			tid=UList[i].peak_idx[k];
			PL[tid].IDX.push_back(curAIdx);
			if(k>=MAX_ISO) PL[tid].IsoTL.push_back(IsoDisT[MAX_ISO-1]);
			else PL[tid].IsoTL.push_back(IsoDisT[k]);
		}		
		curAIdx++;
		if(curAIdx>=pnum) break;

	}	
	ReInitialByLS(ACIs,PL);
	IS=ACIs.size();
	for(i=0;i<IS;i++) ACIs[i].dInt*=ACIs[i].weight;
	return true;
}

void IsoGroup::ReInitialByLS(vector<refPeak> &ACIs,fParList &PL)
{
	//refine the initial coff values
	gsl_matrix *X;
	gsl_vector *c;
	gsl_vector *y;
	size_t curAIdx=ACIs.size();
	if(curAIdx<=0) return;
	size_t pnum=PL.size();
	X=gsl_matrix_alloc(pnum,curAIdx);
	c=gsl_vector_alloc(curAIdx);
	y=gsl_vector_alloc(pnum);
	size_t i,k;
	double *Original;
	Original=new double[curAIdx];
	for(k=0;k<curAIdx;k++)Original[k]=ACIs[k].dInt;
	for(i=0;i<pnum;i++)
	{
		for(k=0;k<curAIdx;k++)
		{
			double x=0;
			int idx=PL[i].FindIdx(k);
			if(idx!=-1) x=PL[i].IsoTL[idx];
			gsl_matrix_set(X,i,k,x);
		}
		gsl_vector_set(y,i,PL[i].Exp_Int);
	}

	if(gsl_robustfit(X,y,c))
	{
		for(i=0;i<curAIdx;i++)
		{
			ACIs[i].dInt=gsl_vector_get(c,i);
			//if(ACIs[i].dInt<0) ACIs[i].dInt=0;
		}
	}
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_vector_free(y);	
	delete []Original;
}

void IsoGroup::GetCharge(vector<int> &ch)
{
	ch.clear();
	if(gPeaks.size()<=0) return;
	size_t sg=gPeaks[0].LocalList.size();
	for(size_t i=0;i<sg;i++) 
	{
		size_t idx=gPeaks[0].LocalList[i].UnitIdx;
		if(UList[idx].charge!=CH_UNKNOWN)ch.push_back(UList[idx].charge);
	}
}
//please make sure the Insnt[0]>eps before fitting
double IsoGroup::Fitting(double gd_cut,int Min_Iso_Num,vector<CalData> &XICPL)
{
	//DecompIP.OutPut();
	size_t iter = 0;
	int status;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;	
	/* Starting point */
	size_t i,k,pnum=gPeaks.size();
	if(pnum<=0) return 0;	
	else if(pnum==1)
	{
		CalData xt;
		xt.ch=CH_UNKNOWN;
		xt.IsoNum=1;
		xt.goodness=1.0;		
		xt.Isomz[0]=gPeaks[0].mz;
		xt.IsotopicE[0]=gPeaks[0].dInt;
		xt.S_N[0]=gPeaks[0].base+gPeaks[0].noise;
		XICPL.push_back(xt);
		return 1.0;
	}
	vector<refPeak> ACIs;
	fParList FTPar;
	//debug
	//if(DecompIP.gPeaks[0].mz>=803.36&&DecompIP.gPeaks[0].mz<803.37)
	//{
	//	int stop=1;
	//	stop++;
	//}
	//
	UpdateWeight(PRE_DEF_WT);//add on 2011.3.29
	SplitCH();

	//printf("  Fitting: initial pars by robust regression:\n" );
	//clock_t start,finish; 
	//start=clock();
	InitialPars(ACIs,FTPar);
	//finish=clock();
	//printf("  End: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);

	//start=finish;
	//printf("  Fitting: begin iteration:\n" );

	gsl_vector *v;		
	size_t vnum=ACIs.size();
	if(vnum<=0) return 0;
	v=gsl_vector_alloc(vnum);
	for(i=0;i<vnum;i++) gsl_vector_set(v,i,ACIs[i].dInt);
	
	gsl_multimin_function_fdf my_func;
	my_func.f = &Res_f;
	my_func.df = &Res_df;
	my_func.fdf = &Res_fdf;
	my_func.n = vnum;
	my_func.params =(void *)(&FTPar);

	T =gsl_multimin_fdfminimizer_conjugate_pr; //gsl_multimin_fdfminimizer_conjugate_fr;//gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc(T,vnum);
	//double stepsize=gsl_vector_get(v,0)/100;
	double stepsize=0.1;
	gsl_multimin_fdfminimizer_set(s,&my_func,v,stepsize,1e-4);
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
		if(status==GSL_ENOPROG)	break;
		status = gsl_multimin_test_gradient(s->gradient,1e-4);
		if (status == GSL_SUCCESS)
		{
			for(i=0;i<vnum;i++)
			{
				double FinialPar=gsl_vector_get(s->x,i);
				ACIs[i].dInt=FinialPar;				
			}
		}
	}while(status==GSL_CONTINUE&&iter<100);
	gsl_multimin_fdfminimizer_free(s);	

	//finish=clock();	
	//printf("  Fitting: End iteration, time=%lf, iter=%d.\n",(double)(finish-start)/CLOCKS_PER_SEC,iter);
	//if fiilure to converage ,using the initial value
	if(iter>=100)
	{
		//printf("\ncatch a failure in fitting!fist m/z=%f, peak num=%d\n",DecompIP.gPeaks[0].mz,DecompIP.gPeaks.size());
		InitialPars(ACIs,FTPar);		
	}

	for(i=0;i<vnum;i++) gsl_vector_set(v,i,ACIs[i].dInt);
	//calulate the goodness of fitting	
	//size_t pnum=FTPar.size();
	double *goodness,*eps;	
	int *IsoNum;	
	goodness=new double[vnum];
	IsoNum=new int[vnum];
	eps=new double[pnum];
	Res_feps(v,(void *)(&FTPar),eps);	
	gsl_vector_free(v);
	for(i=0;i<vnum;i++)
	{		
		if(ACIs[i].dInt<MIN_PEAK_INT) ACIs[i].dInt=0;
		IsoNum[i]=0;	
	}
	for(i=0;i<vnum;i++)
	{
		goodness[i]=0;
		if(ACIs[i].dInt<MIN_PEAK_INT)continue;
		IsoNum[i]=UList[ACIs[i].clusterid].peak_idx.size(); 
		if(IsoNum[i]<=0) continue;	
		double sumi=0;
		for(k=ACIs[i].localid;k<IsoNum[i];k++)
		{
			int tid=UList[ACIs[i].clusterid].peak_idx[k];
			goodness[i]+=fabs(eps[tid]);
			sumi+=gPeaks[tid].dInt;
		}
        IsoNum[i]-=ACIs[i].localid;
		if(IsoNum[i]>MAX_ISO) IsoNum[i]=MAX_ISO;
		if(sumi>1e-2) goodness[i]=1-goodness[i]/sumi;
		else goodness[i]=0.2;
		if(goodness[i]<=0) goodness[i]=0.2;

	/*	if(IsoNum[i]>0) 
		{
			goodness[i]/=IsoNum[i];
			goodness[i]=1-goodness[i];
		}*/
	}	

	//find the max goodness
	double bestG=goodness[0];	
	for(i=1;i<vnum;i++) 
	{
		if(bestG<goodness[i]) bestG=goodness[i];		
	}	
	// generate the calibration data form the iso unit
	CalData xt;	
	for(i=0;i<vnum;i++)
	{
		if(ACIs[i].dInt>MIN_PEAK_INT&&goodness[i]>gd_cut&&IsoNum[i]>=Min_Iso_Num)
		{
			int ui=ACIs[i].clusterid;			
			for(int k=0;k<IsoNum[i];k++)
			{
				xt.Isomz[k]=gPeaks[UList[ui].peak_idx[k]].mz;
				xt.IsotopicE[k]=gPeaks[UList[ui].peak_idx[k]].dInt;
				xt.S_N[k]=gPeaks[UList[ui].peak_idx[k]].base+gPeaks[UList[ui].peak_idx[k]].noise;
				if(xt.S_N[k]>1e-2) xt.S_N[k]=xt.IsotopicE[k]/xt.S_N[k];
			}
			
			xt.ch=UList[ui].charge;
			xt.Isomz[0]=ACIs[i].mz;			
			xt.goodness=goodness[i];
			xt.IsoNum=IsoNum[i];	
			GetIsoDis(xt.Isomz[0]*xt.ch,xt.IsotopicT);
			XICPL.push_back(xt);
		}
	}	
	delete []goodness;
	delete []IsoNum;
	delete []eps;
	return bestG;
}


IsoGroup::IsoGroup()
{

}

IsoGroup::~IsoGroup()
{

}

IsoGroup::IsoGroup(const IsoGroup& cp)
{
	size_t i,sg=cp.gPeaks.size();
	for(i=0;i<sg;i++) 
	{
		gPeaks.push_back(cp.gPeaks[i]);
	}
	sg=cp.UList.size();
	for(i=0;i<sg;i++) 
	{
		UList.push_back(cp.UList[i]);
	}
}

IsoGroup& IsoGroup::operator=(const IsoGroup& cp)
{
	gPeaks.clear();
	size_t i,sg=cp.gPeaks.size();
	for(i=0;i<sg;i++) 
	{
		gPeaks.push_back(cp.gPeaks[i]);
	}
	UList.clear();
	sg=cp.UList.size();
	for(i=0;i<sg;i++) 
	{
		UList.push_back(cp.UList[i]);
	}
	return *this;
}

void IsoGroup::Initial(sPeak &pt)
{	
	gPeaks.clear();
	UList.clear();
	pt.count=1;
	gPeaks.push_back(pt);
	IsoLocal lt;
	lt.UnitIdx=0;
	lt.UnitLocal=0;
	gPeaks[0].LocalList.push_back(lt);
	IsoUnit it;
	it.peak_idx.push_back(0);
	UList.push_back(it);
}

bool IsoGroup::AddPeaks(sPeak &pt,double ppmDM)
{	
	size_t i, gsize=gPeaks.size();
	if(gsize==0)
	{
		Initial(pt);
		return true;
	}

	IsoLocal lt;
	bool IsAdd=false;
	pt.LocalList.clear();

	//printf(">Adding Peaks %d, mz=%lf,ulist=%d\n",gsize,pt.mz, UList.size());

	for(i=0;i<gsize;i++)
	{
		//size_t tt=gPeaks[i].LocalList.size();
		/*printf("Peaks %d\n",i);
			for(size_t kk=0;kk<tt;kk++)
			{
				printf("ch=%d\tpnum=%d\n",UList[gPeaks[i].LocalList[kk].UnitIdx].charge,UList[gPeaks[i].LocalList[kk].UnitIdx].peak_idx.size());
			}*/
		double massdiff=pt.mz-gPeaks[i].mz;
		int ch;
		for(ch=1;ch<=MAX_CH;ch++)
		{
			double isodiff=massdiff-ISO_DIFF[1]/ch;		
			if(fabs(isodiff)/pt.mz<ppmDM/1e6) break;
		}
		if(ch>MAX_CH) continue;
		IsAdd=true;
			
		size_t k, us=gPeaks[i].LocalList.size();
		for(k=0;k<us;k++)
		{					
			if(ch==UList[gPeaks[i].LocalList[k].UnitIdx].charge)//&&gPeaks[i].LocalList[k].UnitLocal==UList[gPeaks[i].LocalList[k].UnitIdx].peak_idx.size()-1) 
			{
				//gPeaks[i].LocalList.push_back(it);
				lt.UnitLocal=UList[gPeaks[i].LocalList[k].UnitIdx].peak_idx.size();
				UList[gPeaks[i].LocalList[k].UnitIdx].peak_idx.push_back(gsize);
				lt.UnitIdx=gPeaks[i].LocalList[k].UnitIdx;	
				pt.LocalList.push_back(lt);
				break;
				
			}
			else if(UList[gPeaks[i].LocalList[k].UnitIdx].charge==CH_UNKNOWN)
			{
				UList[gPeaks[i].LocalList[k].UnitIdx].charge=ch;
				lt.UnitLocal=UList[gPeaks[i].LocalList[k].UnitIdx].peak_idx.size();
				UList[gPeaks[i].LocalList[k].UnitIdx].peak_idx.push_back(gsize);
			    lt.UnitIdx=gPeaks[i].LocalList[k].UnitIdx;	
				pt.LocalList.push_back(lt);
				break;
			}
		}
	    if(k==us)//a new charge
		{
			IsoUnit it;
			it.peak_idx.push_back(i);
			it.peak_idx.push_back(gsize);
			it.charge=ch;
			UList.push_back(it);
			gPeaks[i].count++;//对作为基准峰的次数进行计数
			lt.UnitIdx=UList.size()-1;
			lt.UnitLocal=0;
			gPeaks[i].LocalList.push_back(lt);
			lt.UnitLocal=1;
			pt.LocalList.push_back(lt);
		/*	printf("addtional new, ch=%d\tUid=%d\tPID=%d\t\n",ch,lt.UnitIdx,i);
			printf("\t\tThe peak charge list is:\n");
			size_t tt=gPeaks[i].LocalList.size();
			for(size_t kk=0;kk<tt;kk++)
			{
				printf("ch=%d\tpnum=%d\n",UList[gPeaks[i].LocalList[kk].UnitIdx].charge,UList[gPeaks[i].LocalList[kk].UnitIdx].peak_idx.size());
			}*/
		}	
	}

	if(IsAdd) 
	{		
		gPeaks.push_back(pt);
		//pt.LocalList.clear();
	}
	return IsAdd;
}

//根据峰数目选择，数目越多，排序越靠前，如果数目相同，电荷越低，排序越靠前
void IsoGroup::SortUList()
{
	IsoUnit swp;
	size_t i, sg=UList.size();
	for(i=0;i<sg;i++)
	{
		size_t pnum1=UList[i].peak_idx.size();
		for(size_t k=i+1;k<sg;k++)
		{			
			size_t pnum2=UList[k].peak_idx.size();
			if(pnum1<pnum2)
			{
				swp=UList[i];
				UList[i]=UList[k];
				UList[k]=swp;
			}
			else if(pnum1==pnum2)
			{
				if(UList[i].charge>UList[k].charge)
				{
					swp=UList[i];
					UList[i]=UList[k];
					UList[k]=swp;
				}
			}
		}
	}
}

void IsoGroup::OutPut()
{
	size_t i,sg=UList.size();
	printf("The number of charge cluster is: %d\n",sg);
	for(i=0;i<sg;i++) 
	{
		size_t ts=UList[i].peak_idx.size();
		printf("The peak number for charge= %d is %d\n",UList[i].charge,ts);
		for(size_t k=0;k<ts;k++)
		{
			int tid=UList[i].peak_idx[k];
			printf("\t%lf\t%lf\n",gPeaks[tid].mz,gPeaks[tid].dInt);
		}	
	}
}

double IsoGroup::GetCOS(double *EInt,double *TInt,int num)
{
	int i;
	double Intsum=0;
	for(i=0;i<num;i++) 
	{		
		Intsum+=EInt[i];
	}
	if(Intsum<=0) return 0;
	for(i=0;i<num;i++)EInt[i]=EInt[i]/Intsum;		
	Intsum=0;
	double Intsum1=0;
	double Intsum2=0;
	for(i=0;i<num;i++)
	{
		Intsum+=EInt[i]*EInt[i];
		Intsum1+=TInt[i]*TInt[i];
		Intsum2+=EInt[i]*TInt[i];
	}
	Intsum*=Intsum1;
	if(Intsum<=0) return 0;	
	return Intsum2/sqrt(Intsum);
}


double IsoGroup::GetUnitGD(int ui)
{
	size_t sg=UList.size();
	if(ui>=sg) return 0;
	sg=UList[ui].peak_idx.size();
	if(sg<=0) return 0;
	double *EInt,*TInt;
	EInt=new double[sg];
	TInt=new double[sg];
	for(size_t i=0;i<sg;i++)EInt[i]=gPeaks[UList[ui].peak_idx[i]].dInt;
	GetIsoDis(gPeaks[UList[ui].peak_idx[0]].mz*UList[ui].charge,TInt,sg);
	double breturn=GetCOS(EInt,TInt,sg);
	delete []TInt;
	delete []EInt;
	return breturn;
}


void IsoGroup::UpdateWeight(double Pre_def_w)
{
	size_t i,k,sg=UList.size();
	for(i=0;i<sg;i++) 
	{
		UList[i].weight=1.0;
		for(k=0;k<sg;k++)
		{
			if(k==i) continue;
			if(UList[i].CheckFactor(UList[k]))
			{
				UList[i].weight=Pre_def_w;
				break;
			}
		}
	}
}

void IsoGroup::SplitCH()
{
	size_t k,i,sNum=UList.size();
	if(sNum<=0) return;
	size_t pnum=gPeaks.size();
	if(pnum<=0) return;

	int tid;
	int curAIdx=0;
	size_t IS;
	double IsoDisT[MAX_ISO];
	IsoUnit IT;
	IT.weight=1.0;
	//add the charge determined overlap at first
	for(i=0;i<sNum;i++)
	{	
		if(UList[i].weight<0.9) continue;
		IS=UList[i].peak_idx.size();
		if(IS<=0) continue;
		tid=UList[i].peak_idx[0];
		double mass=gPeaks[tid].mz*UList[i].charge;
		GetIsoDis(mass,IsoDisT);

		for(k=1;k<IS;k++)
		{
			int tidnew=UList[i].peak_idx[k];		
			if(IsANew(mass,gPeaks[tid].dInt,gPeaks[tidnew].dInt,k-1))
			{	
				if(k==IS-1) IT.charge=CH_UNKNOWN;//the last one can not determine the charge
				else IT.charge=UList[i].charge;
				gPeaks[tidnew].count++;//基准峰计数

				for(size_t j=k;j<IS;j++)
				{
					if(GetIsoDis(mass,j-k+1)<1e-4) break;
					int tidx=UList[i].peak_idx[j];	
					IT.peak_idx.push_back(tidx);					
				}	
				UList.push_back(IT);//注意，新加入的UI不会被分裂算法考虑
				IT.peak_idx.clear();				
			}
			tid=tidnew;
		}
		UList[i].AdjustCH(mass);
	}
}

bool IsoGroup::IsANew(double mass, double Iso1,double Iso2,int idx)
{
	double LR;
	if(idx==0)LR=0.000584*mass+0.2015;
	else if(idx==1) LR=0.0003072*mass+0.1738;
	else if(idx==2) LR=0.0001991*mass+0.2665;
	else if(idx==3) LR=0.0001342*mass+0.4034;
	else  LR=0.0001557*mass+0.1536;//all >6 is the same with 6
	if(Iso2>LR*Iso1*FACTOR_MUL) return true;
	return false;
}


//---------------------------------------------------------------------------//
//--------------default restrictions without check---------------------------//
double Res_f(const gsl_vector *v, void *params)
{
	fParList *Par;
	Par=(fParList *)params;
	size_t i,IntsNum=Par->size();
	if(IntsNum<=0) return NAN; //if no data return the maximal float
	double *eps;
	eps=new double[IntsNum];		
	Res_feps(v, params,eps);
	double sum=0;	
	for(i=0;i<IntsNum;i++) 
	{
		sum+=eps[i]*eps[i];	
	}

	size_t par_num=v->size;
	double total_Int=0;
	for(i=0;i<par_num;i++)
	{
		double tmpf=gsl_vector_get(v,i);		
		// for the negative values of C
		if(tmpf<-ABS_SIG_CUT) sum+=F_K*tmpf*tmpf;
		//add on 2009.4.25
		//for to many non-zero parameters		
		//else if(tmpf>ABS_SIG_CUT) total_Int+=tmpf;//tmpf*tmpf;
		//added on 2011.1.11
		else if(tmpf>ABS_SIG_CUT) total_Int++;
	}		
	sum+=total_Int*SIGNAL_FACTOR;
	delete []eps;
	return sum;
}

void Res_feps(const gsl_vector *v, void *params,double *eps)
{
	fParList *Par;
	Par=(fParList *)params;
	size_t IntsNum=Par->size();
	size_t i,k;	
	for(i=0;i<IntsNum;i++) eps[i]=Par->at(i).Exp_Int;		
	for(i=0;i<IntsNum;i++)
	{
		size_t sNum=Par->at(i).IsoTL.size();
		for(k=0;k<sNum;k++)
		{
			eps[i]-=Par->at(i).IsoTL[k]*gsl_vector_get(v,Par->at(i).IDX[k]);
		}
	}
}

//--------------default restrictions without check---------------------------//
void Res_df(const gsl_vector *v, void *params,gsl_vector *df)
{	
	fParList *Par;
	Par=(fParList *)params;
	size_t i,k,IntsNum=Par->size();
	double *eps;
	eps=new double[IntsNum];
	Res_feps(v,params,eps);
	size_t par_num=v->size;	
	for(i=0;i<par_num;i++)
	{
		double dfP=0;
		for(k=0;k<IntsNum;k++)
		{
			int idx=Par->at(k).FindIdx(i);
			if(idx!=-1) 
			{
				dfP+=-2*eps[k]*Par->at(k).IsoTL[idx];
			}
		}		
		//add on 2009.4.25 for the negative values of C
		double tmpf=gsl_vector_get(v,i);
		if(tmpf<-ABS_SIG_CUT) dfP+=2*F_K*tmpf;
		//else if(tmpf>ABS_SIG_CUT) dfP+=1;//dfP+=2*tmpf;//modified on 2011.5.7
		////////////////////////////////////////
		gsl_vector_set(df,i,dfP);
	}
	delete []eps;
}	

void Res_fdf(const gsl_vector *v, void *params,double *f,gsl_vector *df)
{	
	*f=Res_f(v,params);
	Res_df(v,params,df);
}
