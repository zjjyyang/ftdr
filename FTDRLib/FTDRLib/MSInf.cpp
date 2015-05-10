#include "StdAfx.h"
#include "MSInf.h"
#include "GaussFit.h"

MSInf::MSInf(void)
{
	pkl=NULL;
	PNum=0;
	isCentroided=false;
	Signal_Th=1e-4;
	RSlower=0.8;
	Min_Iso_Num=2;
	GD_CUT=0.01;
	MassErrorTol=10;
}

MSInf::~MSInf(void)
{
	if(PNum>0&&pkl!=NULL) delete []pkl;
}

MSInf::MSInf(const MSInf &r)
{
	isCentroided=r.isCentroided;
	Signal_Th=r.Signal_Th;
	RSlower=r.RSlower;
	PNum=r.PNum;
	ScanNumber=r.ScanNumber;
	RT=r.RT;
	if(PNum<=0) return;
	pkl=new DataPeak[PNum];	
	memcpy((void *)pkl,(void *)(r.pkl),sizeof(double)*4*PNum);
}

MSInf & MSInf::operator =(const MSInf &r)
{
	if(PNum>0&&pkl!=NULL) 
	{
		delete []pkl;
		pkl=NULL;
		PNum=0;
	}
	isCentroided=r.isCentroided;
	Signal_Th=r.Signal_Th;
	RSlower=r.RSlower;
	ScanNumber=r.ScanNumber;
	RT=r.RT;
	PNum=r.PNum;
	if(PNum<=0) return *this;
	pkl=new DataPeak[PNum];
	memcpy((void *)pkl,(void *)(r.pkl),sizeof(double)*4*PNum);
	return *this;
}

bool MSInf::SetPKL(double *pdval,int dim)
{
	if(PNum>0&&pkl!=NULL) 
	{
		delete []pkl;
		pkl=NULL;
		PNum=0;
	}
	PNum=dim;
	if(PNum<=0) return false;
	pkl=new DataPeak[PNum];
	for (int inx = 0; inx < dim; inx++)
	{
		pkl[inx].dMass=  (double)pdval[((inx)*6)+0] ;
		pkl[inx].dIntensity	=  (double)	pdval[((inx)*6)+1] ;
		//fRes		=  (float) 	pdval[((inx)*6)+2] ;
		pkl[inx].dBase		=  (float) 	pdval[((inx)*6)+3] ;
		pkl[inx].dNoise		=  (float) 	pdval[((inx)*6)+4] ;		
	}	
	return true;
}

bool MSInf::SetPKL(vector<double> &mz,vector<double> &abu)
{
	isCentroided =false;
	if(PNum>0&&pkl!=NULL) 
	{
		delete []pkl;
		pkl=NULL;
		PNum=0;
	}
	PNum=mz.size();
	int count=abu.size();
	if(PNum>count) PNum=count;
	if(PNum<=0) return false;
	pkl=new DataPeak[PNum];
	for(size_t i=0;i<PNum;i++)
	{
		pkl[i].dMass=mz[i];
		pkl[i].dIntensity=abu[i];
	}
	return true;
}

bool MSInf::SetPKL(DataPeak *tmpPKL,int count)
{
	if(PNum>0&&pkl!=NULL) 
	{
		delete []pkl;
		pkl=NULL;
		PNum=0;
	}
	PNum=count;
	if(PNum<=0) return false;
	pkl=new DataPeak[PNum];
	memcpy((void *)pkl,(void *)(tmpPKL),sizeof(double)*2*PNum);
	return true;
}

int MSInf::Centriod(string instrument)
{
	if(PNum<=0) return 0;
	if(isCentroided) return 0;

	double dm;
	int begin=0;	
	vector<double> c_intarray;
	vector<double> c_mzarray;
	int i;
	while(begin<PNum)
	{				
		int end=begin+1;
		int idx_old=begin;
		dm=1.5*GetSampleDm(instrument,pkl[begin].dMass);
		while(end<PNum)
		{
			if(pkl[end].dIntensity<Signal_Th) break;
			if(pkl[end].dMass>pkl[idx_old].dMass+dm) break;			
			idx_old=end;
			end++;
		}
		vector<double> pkmass;
		vector<double> pkintensity;
		vector<double> beta;
		for(i=begin;i<end;i++)
		{	
			pkmass.push_back(pkl[i].dMass);
			pkintensity.push_back(pkl[i].dIntensity);
		}
		double RSquare=Gauss_Fit(pkmass,pkintensity,beta);
		if(RSquare>RSlower)
		{
			int n=beta.size()/4;
			int mu_idx;
			int abu_idx;
			for(i=0;i<n;i++)
			{
				abu_idx=4*i;
				mu_idx=abu_idx+1;
				if(beta[abu_idx]<=0||beta[mu_idx]<=0) continue;
				c_mzarray.push_back(beta[mu_idx]);
				c_intarray.push_back(beta[abu_idx]);
			}
		}
		begin=end;
	}

	delete []pkl;
	pkl=NULL;
	PNum=c_intarray.size();
	if(PNum<=0) return 0;
	pkl=new DataPeak[PNum];
	for(i=0;i<PNum;i++)
	{
		pkl[i].dIntensity=c_intarray[i];
		pkl[i].dMass=c_mzarray[i];
	}	
	isCentroided= true;
	return PNum;
}

bool MSInf::OutPutData(FILE *fp)
{
	if(fp==NULL) return false;
	fprintf(fp,">%d.%d\n",ScanNumber,PNum);
	for(int i=0;i<PNum;i++)
	{
		fprintf(fp,"%lf\t%lf\n",pkl[i].dMass,pkl[i].dIntensity);
	}
	return true;
}

double MSInf::GetSampleDm(string instrument,double mass)
{
	if(instrument=="FT") return 6.06e-9*mass*mass-7.063e-11*mass+4.085e-8;//for FT
	else if(instrument=="Orbitrap") return 1.12e-009*mass*mass+2.204e-006*mass-0.0003366;//for Orbitrap
	else return 0;
}


void MSInf::RemUsed(vector<Peak> &tmpPKL)
{
	size_t i,pnum=tmpPKL.size();
	vector<Peak> _tmppkl;
	for(i=0;i<pnum;i++)
	{
		if(tmpPKL[i].IsUsed) continue;
		_tmppkl.push_back(tmpPKL[i]);
	}
	
	pnum=_tmppkl.size();
	tmpPKL.clear();
	for(i=0;i<pnum;i++)
	{
		tmpPKL.push_back(_tmppkl[i]);
	}
}

int MSInf::FindIsoCPmz(double BaseInt,FILE *fp,vector<Peak> &tmpPKL)
{
	size_t pnum=tmpPKL.size();
	if(pnum<Min_Iso_Num) return NO_PEAKS;
	double oldmz=tmpPKL[0].mz;
	tmpPKL[0].IsUsed=true;
	if(tmpPKL[0].dInt<RINT_CUT*BaseInt)
	{	
		RemUsed(tmpPKL);		
		return CPEAK_FAILUR;
	}
	vector<XICPeak> tmpXICPL;
	IsoGroup IG;
	size_t i;	
	IG.Initial(tmpPKL[0]);
	for(i=1;i<pnum;i++) 
	{
		if(IG.AddPeaks(tmpPKL[i],MassErrorTol))
		{
			tmpPKL[i].IsUsed=true;
			oldmz=tmpPKL[i].mz;
		}	
		else if(tmpPKL[i].mz>oldmz+1.1) break;			
	}

	if(IG.UList.size()>1)
	{
		printf("\nBefore validation:find a mutiple charge overlap case:scannum=%d\tbasemz=%lf\n",ScanNumber,IG.gPeaks[0].mz);
	}
	
	if(IG.gPeaks.size()<Min_Iso_Num)//not a iso group, remove these peaks only
	{
		RemUsed(tmpPKL);		
		return CPEAK_FAILUR;
	}
	double CurrentBestGD=Fitting(IG,GD_CUT,Min_Iso_Num,tmpXICPL);
	if(CurrentBestGD<GD_CUT) 
	{
		RemUsed(tmpPKL);
		return NO_GOOD_CLUSTER;	
	}
	size_t si=tmpXICPL.size();	
	vector<int> chList;
	for(i=0;i<si;i++)
	{
		if(BaseInt>DEPSL) tmpXICPL[i].RelInt=tmpXICPL[i].MonoInt/BaseInt;
		tmpXICPL[i].ScanNum=ScanNumber;		
		tmpXICPL[i].Out2FileNew(fp);
		size_t k,chNum=chList.size();
		for(k=0;k<chNum;k++)
		{
			if(chList[k]==tmpXICPL[i].charge) break;
		}
		if(k==chNum) chList.push_back(tmpXICPL[i].charge);
	}

	if(chList.size()>1)
	{
		printf("\After validation: a mutiple charge case:charger_num=%d\tscannum=%d\tbasemz=%lf\n",chList.size(),ScanNumber,IG.gPeaks[0].mz);
	}

	RemUsed(tmpPKL);
	return FIND_ONE;
}

double MSInf::GetBasePeak()
{
	double MaxInt=0;
	for(int i=0;i<PNum;i++)
	{
		if(MaxInt<pkl[i].dIntensity) MaxInt=pkl[i].dIntensity;
	}
	return MaxInt;
}

int MSInf::FindByPmz(FILE *fp)
{
	if(PNum<=0) return 0;
	vector<Peak> tmpPKL;	
	Peak pt;
	double baseInt=GetBasePeak();
	for(int i=0;i<PNum;i++) 
	{
		pt=pkl[i];
		tmpPKL.push_back(pt);
	}
	int CT=0;	
	int bT=FindIsoCPmz(baseInt,fp,tmpPKL);
	if(bT==FIND_ONE) CT++;
	while(bT!=NO_PEAKS) 
	{	
		bT=FindIsoCPmz(baseInt,fp,tmpPKL);
		if(bT==FIND_ONE)CT++;		
	}	
	return CT;
}