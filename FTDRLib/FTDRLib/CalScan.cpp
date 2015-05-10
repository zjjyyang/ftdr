#include "stdafx.h"
#include "CalScan.h"
#include "math.h"


CalScan::CalScan(void)
{	
	ISt=0;
	ESt=0;
	ScanNum=0;
	RT=0;
	Iso_win=2.0;//default isolation window is 2.0Da
	Pre_dm=10.0;//ppm
	Min_Iso_Num=1;
	BaseInt=0;
}


CalScan::~CalScan(void)
{

}

CalScan::CalScan(const CalScan &cp)
{	
	size_t i,sg=cp.PKL.size();	
	for(i=0;i<sg;i++)
	{
		PKL.push_back(cp.PKL[i]);
	}
	sg=cp.StatusPar.size();
	for(i=0;i<sg;i++) StatusPar.push_back(cp.StatusPar[i]);
	ISt=cp.ISt;
	ESt=cp.ESt;
	ScanNum=cp.ScanNum;
	RT=cp.RT; 
	Iso_win=cp.Iso_win;
	Pre_dm=cp.Pre_dm;
	BaseInt=cp.BaseInt;
	Min_Iso_Num=cp.Min_Iso_Num;
}


CalScan& CalScan::operator=(const CalScan &cp)
{
	PKL.clear();
	size_t i,sg=cp.PKL.size();	
	for(i=0;i<sg;i++)
	{
		PKL.push_back(cp.PKL[i]);
	}

	sg=cp.StatusPar.size();
	StatusPar.clear();
	for(i=0;i<sg;i++) StatusPar.push_back(cp.StatusPar[i]);
	ISt=cp.ISt;
	ESt=cp.ESt;
	ScanNum=cp.ScanNum;
	RT=cp.RT;
	Iso_win=cp.Iso_win;
	Pre_dm=cp.Pre_dm;
	BaseInt=cp.BaseInt;
	Min_Iso_Num=cp.Min_Iso_Num;
	return *this;
}

bool CalScan::Convert(Scan *cs,vector<size_t> &PreSuvList)
{
	if(cs==NULL) return false;
	RT=cs->retentionTimeInSec_/60;
	ScanNum=cs->curScanNum;
	size_t i, TS=cs->status_par_name.size();
	size_t ss=PreSuvList.size();
	StatusPar.clear();
	bool IsSecond=false;
	for(i=0;i<ss;i++)
	{
		size_t tid=PreSuvList[i];
		StatusPar.push_back(atof(cs->status_par_value[tid].c_str()));
	}

	ss=cs->tailer_par_name.size();
	for(i=0;i<ss;i++)
	{
		if(cs->tailer_par_name[i].find("Ion Injection Time")!=string::npos)
		{
			ISt=atof(cs->tailer_par_value[i].c_str());
		}
		else if(cs->tailer_par_name[i].find("Elapsed Scan Time")!=string::npos)
		{
			ESt=atof(cs->tailer_par_value[i].c_str());
		}
	}
	PKL.clear();
	ss=cs->getNumDataPoints();
	mPeak pt;
	for(i=0;i<ss;i++)
	{
		pt.dMass=cs->mzArray_[i];
		pt.dIntensity=cs->intensityArray_[i];
		pt.rInt=pt.dIntensity/cs->basePeakIntensity_;
		pt.noise=cs->noise_[i];
		pt.baseLine=cs->baseline_[i];
		PKL.push_back(pt);
	}
	BaseInt=cs->basePeakIntensity_;
	return true;
}


bool CalScan::Convert(Scan *cs,vector<size_t> &PreSuvList,size_t TailerIDX[2])
{
	if(cs==NULL) return false;
	RT=cs->retentionTimeInSec_/60;
	ScanNum=cs->curScanNum;
	size_t i;// TS=cs->status_par_name.size();
	size_t ss=PreSuvList.size();
	//StatusPar.resize();
	StatusPar.clear();
	for(i=0;i<ss;i++)
	{
		size_t tid=PreSuvList[i];
		StatusPar.push_back(atof(cs->status_par_value[tid].c_str()));
	}

	ISt=atof(cs->tailer_par_value[TailerIDX[0]].c_str());
	ESt=atof(cs->tailer_par_value[TailerIDX[1]].c_str());	
	
	ss=cs->getNumDataPoints();
	//mPeak pt;
	PKL.resize(ss);		
	for(i=0;i<ss;i++)
	{
		PKL[i].dMass=cs->mzArray_[i];
		PKL[i].dIntensity=cs->intensityArray_[i];
		PKL[i].rInt=PKL[i].dIntensity/cs->basePeakIntensity_;
		PKL[i].noise=cs->noise_[i];
		PKL[i].baseLine=cs->baseline_[i];		
	}
	BaseInt=cs->basePeakIntensity_;
	return true;
}

void CalScan::InitialSTSisze(size_t sg)
{
	StatusPar.resize(sg);
}

//you must initial StatusPar by function InitialStatus;
//bool CalScan::Convert(Scan *cs,vector<size_t> &PreSuvList,size_t TailerIDX[2])
//{
//	if(cs==NULL) return false;
//	RT=cs->retentionTimeInSec_/60;
//	ScanNum=cs->curScanNum;
//	size_t i;
//	size_t ss=PreSuvList.size();
//	for(i=0;i<ss;i++)
//	{
//		size_t tid=PreSuvList[i];
//		StatusPar[i]=(atof(cs->status_par_value[tid].c_str()));
//	}
//
//	ISt=atof(cs->tailer_par_value[TailerIDX[0]].c_str());
//	ESt=atof(cs->tailer_par_value[TailerIDX[1]].c_str());	
//	
//	ss=cs->getNumDataPoints();
//	mPeak pt;
//	PKL.resize(ss);		
//	for(i=0;i<ss;i++)
//	{
//		PKL[i].dMass=cs->mzArray_[i];
//		PKL[i].dIntensity=cs->intensityArray_[i];
//		PKL[i].rInt=pt.dIntensity/cs->basePeakIntensity_;
//		PKL[i].noise=cs->noise_[i];
//		PKL[i].baseLine=cs->baseline_[i];		
//	}
//	BaseInt=cs->basePeakIntensity_;
//	return true;
//}


//the user must make sure the index no t exceed the band
int CalScan::TFindL(int b,int e,double fmz)
{
	if(e<b) return -1;
	double dm=Pre_dm*fmz/1e6;
	if(PKL[e].dMass<fmz-dm) return -1;
	if(PKL[b].dMass>fmz+dm) return -1;
	if(e-b<=1)
	{
		//if(PKL[b].dMass>fmz+dm||PKL[b].dMass<fmz-dm) return -1;
		return b;	
	}
	int idx=(b+e)/2;	
	if(PKL[idx].dMass>fmz)
	{
		return TFindL(b,idx,fmz);
	}
	else return TFindL(idx,e,fmz);
}

int CalScan::TFindH(int b,int e,double fmz)
{
	if(e<b) return -1;
	double dm=Pre_dm*fmz/1e6;
	if(PKL[e].dMass<fmz-dm) return -1;
	if(PKL[b].dMass>fmz+dm) return -1;
	if(e-b<=1)
	{
		//if(PKL[e].dMass>fmz+dm||PKL[e].dMass<fmz-dm) return -1;
		return e;
	}
	int idx=(b+e)/2;
	if(PKL[idx].dMass>fmz)
	{
		return TFindH(b,idx,fmz);
	}
	else return TFindH(idx,e,fmz);
}

//Error less ciration
//lasted modifide on 2012.2.7
int CalScan::LocatePmz(double DM,double pmz)
{
	int middle=-1;
	int n=PKL.size();
	if(n<=0) return -1;
	int iB=TFindL(0,n-1,pmz);
	if(iB==-1) return -1;
	double dm1,dm2;
	dm1=fabs(PKL[iB].dMass-pmz);
	if(iB<n-1) dm2=fabs(PKL[iB+1].dMass-pmz);
	if(dm1>dm2)
	{
		middle=iB+1;
		dm1=dm2;
	}
	else middle=iB;
	if(dm1>DM) return -1;	
	return middle;
}

//Get the most high in the limit range
//lasted modifide on 2012.2.7,retrun to the original
int CalScan::LocatePmz(double begin,double end,double pmz)
{
	int i,middle=-1;	
	int n=PKL.size();
	if(n<=0) return -1;
	double MaxInt=0;

	int iB=TFindL(0,n-1,begin);
	if(iB==-1) return -1;
	int iE=TFindH(0,n-1,end);
	if(iE==-1) return -1;
	//double dm=Pre_dm*pmz/1e6;

	for(i=iB;i<=iE;i++)
	{		
		if(PKL[i].dIntensity<=1e-6) continue;	
	/*	if(pmz>PKL[i].dMass+dm) continue;
		if(pmz<PKL[i].dMass-dm) break;*/
		if(PKL[i].dMass<begin) continue;
		if(PKL[i].dMass>end) break;
		if(MaxInt<PKL[i].dIntensity)
		{
			MaxInt=PKL[i].dIntensity;
			middle=i;
		}
	}
	return middle;
}

//get a single peak by mz values, find the most high one in a dm range
bool CalScan::Findmz(mPeak &mt,double pmz,double dm)
{	
	int middle=LocatePmz(pmz-dm,pmz+dm,pmz);
	if(middle==-1) return false; //can not find, use the original	
	mt=PKL[middle];	
	return true;
}

double CalScan::GetTIC()
{
	size_t sg=PKL.size();
	double tic=0;
	for(size_t i=0;i<sg;i++) tic+=PKL[i].dIntensity/TIC_SCALE;
	return tic;
}

void CalScan::PackageIt(CalData &it)
{
	it.BaseInt=BaseInt;
	it.RT=RT;
	//it.IIt=ISt;
	//it.ESt=ESt;
	it.scannum=ScanNum;
	it.tic=GetTIC();

	it.FTStatus.clear();
	it.FTStatus.push_back(ISt);
	it.FTStatus.push_back(ESt);
	size_t sg=StatusPar.size();
	for(int i=0;i<sg;i++)
	{
		it.FTStatus.push_back(StatusPar[i]);
	}	

}
//Extract  all ISO peaks by pmz, it is the mono-mz
//do not take into account the pmass lost case
//that case was processed in the first finding as a special case
bool CalScan::FindISO(CalData &it,int ch,double pmz)
{	
	return FindISO_dMET(it,ch,pmz,-Pre_dm,Pre_dm);
}

//Extract  all ISO peaks by pmz, it is the mono-mz
//do not take into account the pmass lost case
//that case was processed in the first finding as a special case
bool CalScan::FindISO_dMET(CalData &it,int ch,double pmz,double ppmL,double ppmH)
{	
	if(ch<0) return false;	
	//if(IsoDisT[0]<=1e-4) return false;

	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		it.IsotopicE[i]=0;
		it.Isomz[i]=0;	
		it.S_N[i]=0;
	}

	PackageIt(it);	
	it.ch=ch;
	GetIsoDis(pmz*ch,it.IsotopicT);	
	double fmz=pmz;	
	//a bug find on 2012.2.19
	double dmL=fmz*ppmL/1e6;
	double dmH=fmz*ppmH/1e6;
	///
	double sumI=0;
	for(i=0;i<MAX_ISO;i++) sumI+=it.IsotopicT[i];

	if(ch==CH_UNKNOWN)
	{
		int idx=LocatePmz(fmz+dmL,fmz+dmH,fmz);
		if(idx==-1) return false;
		it.S_N[0]=PKL[idx].noise+PKL[idx].baseLine;
		it.Isomz[0]=PKL[idx].dMass;
		it.IsotopicE[0]=PKL[idx].dIntensity;	
		it.IsoNum=1;
		if(it.S_N[0]>1e-2)it.S_N[0]=it.IsotopicE[0]/it.S_N[0];
		else it.S_N[0]=0;
		if(it.IsotopicE[0]<RINT_CUT*BaseInt) return false;	
		if(it.S_N[0]<MIN_SN) return false;
		it.goodness=it.IsotopicT[0]/sumI;
	}
	else
	{
		double sumE=0;
		double sumIE=0;
		double gd,gdold=0;
		for(i=0;i<MAX_ISO;i++)
		{
			fmz=pmz+ISO_DIFF[i]/ch;
			int idx=LocatePmz(fmz+dmL,fmz+dmH,fmz);
			if(idx==-1) break;
			it.S_N[i]=PKL[idx].noise+PKL[idx].baseLine;
			it.Isomz[i]=PKL[idx].dMass;
			it.IsotopicE[i]=PKL[idx].dIntensity;
			if(it.S_N[i]>1e-2)it.S_N[i]=it.IsotopicE[i]/it.S_N[i];
			else it.S_N[i]=0;
			if(it.IsotopicE[i]<RINT_CUT*BaseInt)break;	
			if(it.S_N[i]<MIN_SN) break;			
			sumE+=it.IsotopicE[i];
			sumIE+=it.IsotopicT[i]*it.IsotopicE[i];				
			if(it.IsotopicT[i]<0.5*it.IsotopicT[0]) break;
			double ratioT=it.IsotopicT[i]/it.IsotopicT[0];
			double ratioE=it.IsotopicE[i]/it.IsotopicE[0];
			if(ratioE>1.5*ratioT||ratioE<0.5*ratioT) break;
		}		
		it.IsoNum=i;
		if(i<ISO_CUT) return false;	
		it.goodness=sumIE/(sumI*sumE);
	}	
	if(it.goodness<GD_CUT) return false;
	return true;
}

void CalScan::RemUsed(vector<sPeak> &tmpPKL)
{
	size_t i,pnum=tmpPKL.size();
	vector<sPeak> _tmppkl;
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

int CalScan::FindIsoCPmz(vector<sPeak> &tmpPKL,vector<CalData> &PIso)
{	
	size_t i, pnum=tmpPKL.size();
	if(pnum<Min_Iso_Num) return NO_PEAKS;
	double oldmz=tmpPKL[0].mz;
	tmpPKL[0].IsUsed=true;		
	IsoGroup IG;
	IG.Initial(tmpPKL[0]);
	for(i=1;i<pnum;i++) 
	{
		if(IG.AddPeaks(tmpPKL[i],Pre_dm))
		{
			tmpPKL[i].IsUsed=true;
			oldmz=tmpPKL[i].mz;
		}	
		else if(tmpPKL[i].mz>oldmz+1.1) break;	
		//printf("the unist validate: %d\n",IG.UList.size());
	}
	
	if(IG.gPeaks.size()<Min_Iso_Num)
	{
		RemUsed(tmpPKL);
		return CPEAK_FAILUR;
	}

	double CurrentBestGD=IG.Fitting(GD_CUT_L,Min_Iso_Num,PIso);//GD_CUT
	if(CurrentBestGD<GD_CUT_H) 
	{
		RemUsed(tmpPKL);
		return NO_GOOD_CLUSTER;	
	}
	//size_t si=PIso.size();	
	//for(i=0;i<si;i++)PIso[i].BaseInt=BaseInt;	
	RemUsed(tmpPKL);
	return FIND_ONE;
}

//try to find the mono-pmz from a initial mz and ch, and try to find all possible pmzs 
//in the isolotion window, generally it is +/-2.0Da
//The IW will be extracted from the raw file directly
bool CalScan::ExtendFindISO(vector<CalData> &PIso,double pmz,double IW)
{
	//for debug
	//if(ScanNum==5919)
	//{
	//	ScanNum=5919;
	//}
	//end
	double CurrentBestGD=0;
	vector<sPeak> tmpPKL;
	size_t pnum=PKL.size();
	if(pnum<=0) return false;
	int begin=TFindL(0,pnum-1,pmz-IW);
	if(begin==-1) return false;
	while(begin<pnum)
	{
		if(PKL[begin].dMass>=pmz-IW) break;
		begin++;
	}
	int end=TFindH(0,pnum-1,pmz+IW);
	if(end==-1) return false;
	while(end>begin)
	{
		if(PKL[end].dMass<=pmz+IW) break;
		end--;
	}

	sPeak pt;
	vector<sPeak> tmpPKL1;

	int k,i=begin-1;	
	//extend according the iso pair rules
	bool IsAdd=false;
	double K=Pre_dm/1e6;
	for(;i>0;i--)
	{
		if(PKL[begin].dMass-PKL[i].dMass>ISO_DIFF[1]) break;
		for(k=begin;k<=end;k++)
		{
			double delt=PKL[k].dMass-PKL[i].dMass;
			int ch;
			for(ch=MAX_CHARGE;ch>=1;ch--)
			{
				if(fabs(delt*ch-ISO_DIFF[1])<K*PKL[i].dMass)
				{
					IsAdd=true;
					break;
				}
			}
		}
		if(IsAdd) 
		{
			//begin=i;//adjust the begin
			pt=PKL[i];
			tmpPKL1.push_back(pt);
			IsAdd=false;
		}
	}
	k=tmpPKL1.size()-1;
	for(;k>=0;k--)
	{
		tmpPKL.push_back(tmpPKL1[k]);
	}
	for(i=begin;i<=end;i++) 
	{
		pt=PKL[i];
		tmpPKL.push_back(pt);
	}
	i=end+1;
	for(;i<pnum;i++)
	{
		if(PKL[i].dMass-PKL[end].dMass>ISO_DIFF[1]) break;
		for(k=begin;k<=end;k++)
		{
			double delt=PKL[i].dMass-PKL[k].dMass;
			int ch;
			for(ch=MAX_CHARGE;ch>=1;ch--)
			{
				if(fabs(delt*ch-ISO_DIFF[1])<K*PKL[i].dMass)
				{
					IsAdd=true;
					break;
				}
			}
		}
		if(IsAdd) 
		{
			//end=i;//adjust the begin
			pt=PKL[i];
			tmpPKL.push_back(pt);
			IsAdd=false;
		}
	}
	
	//int CT=0;	
	int bT=FindIsoCPmz(tmpPKL,PIso);
	//if(bT==FIND_ONE)CT++;
	while(bT!=NO_PEAKS) 
	{	
		bT=FindIsoCPmz(tmpPKL,PIso);
		//if(bT==FIND_ONE)CT++;
		//printf("\r%d",CT);
	}	
	pnum=PIso.size();	
	//here more filteration may be provide for t he addtional pmzs
	vector<CalData> tmpPIso;
	for(i=0;i<pnum;i++)
	{
		if(!PIso[i].IspmzIn(pmz,Pre_dm))
		{
			if(PIso[i].GetIsoNum()>MIN_H_ISO_NUM
				&&PIso[i].IsotopicE[0]>BaseInt*R_INT_CUT_H
				&&PIso[i].ch!=CH_UNKNOWN
				&&PIso[i].goodness>GD_CUT_H) 
			{
				tmpPIso.push_back(PIso[i]);
			}

		}
		else tmpPIso.push_back(PIso[i]);
	}
	/////////////////
	PIso.clear();
	pnum=tmpPIso.size();
	for(i=0;i<pnum;i++) PIso.push_back(tmpPIso[i]);
	for(i=0;i<pnum;i++)	PackageIt(PIso[i]);
	return true;	
}

void CalScan::RetrivalPeaks(double min_mz,double max_mz,vector<mPeak> &tmpPKL)
{	
	int n=PKL.size()-1;
	int bg=TFindL(0,n,min_mz);
	if(bg==-1) return;
	int ed=TFindL(0,n,max_mz);
	if(ed==-1) return;
	for(int i=bg;i<=ed;i++) tmpPKL.push_back(PKL[i]);		
}

//whole spectrum finding
bool CalScan::FindALLISO(vector<CalData> &PIso)
{
	double CurrentBestGD=0;
	vector<sPeak> tmpPKL;
	size_t i,pnum=PKL.size();
	if(pnum<=0) return false;	
	
	sPeak pt;
	vector<sPeak> tmpPKL1;	
	for(i=0;i<pnum;i++)
	{
		pt=PKL[i];
		tmpPKL.push_back(pt);
	}

	int bT=FindIsoCPmz(tmpPKL,PIso);	
	while(bT!=NO_PEAKS) 
	{	
		bT=FindIsoCPmz(tmpPKL,PIso);
	}	
	pnum=PIso.size();	
	//here more filteration may be provide for t he addtional pmzs	
	for(i=0;i<pnum;i++)	PackageIt(PIso[i]);
	return true;	
}

//will return the file pos of next spectrum
fpos_t CalScan::Write2File(FILE *fp)
{	
	//write fixed terms
	fwrite(&ScanNum,sizeof(long),1,fp);
	fwrite(&RT,sizeof(double),1,fp);
	fwrite(&ISt,sizeof(double),1,fp);
	fwrite(&ESt,sizeof(double),1,fp);
	fwrite(&BaseInt,sizeof(double),1,fp);
	fwrite(&Iso_win,sizeof(double),1,fp);
	size_t i,sg=StatusPar.size();
	fwrite(&sg,sizeof(size_t),1,fp);
	for(i=0;i<sg;i++)
	{
		fwrite(&(StatusPar[i]),sizeof(double),1,fp);

	}
	sg=PKL.size();
	fwrite(&sg,sizeof(size_t),1,fp);
	for(i=0;i<sg;i++)
	{
		PKL[i].Write2File(fp);
	}
	fpos_t fposB;
	fgetpos(fp,&(fposB));
	return fposB;
}

//will return the file pos of next spectrum
void CalScan::ReadFromFile(FILE *fp)
{	
	//write fixed terms
	fread(&ScanNum,sizeof(long),1,fp);
	fread(&RT,sizeof(double),1,fp);
	fread(&ISt,sizeof(double),1,fp);
	fread(&ESt,sizeof(double),1,fp);
	fread(&BaseInt,sizeof(double),1,fp);
	fread(&Iso_win,sizeof(double),1,fp);
	size_t i,sg;
	StatusPar.clear();
	fread(&sg,sizeof(size_t),1,fp);
	double tmpf;
	for(i=0;i<sg;i++)
	{
		fread(&tmpf,sizeof(double),1,fp);
		StatusPar.push_back(tmpf);
	}
	fread(&sg,sizeof(size_t),1,fp);
	PKL.clear();
	mPeak pt;
	for(i=0;i<sg;i++)
	{
		pt.ReadFromFile(fp);
		PKL.push_back(pt);
	}	
}