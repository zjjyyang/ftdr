#include "StdAfx.h"
#include "IdxRaw.h"
#include "math.h"


IdxRaw::IdxRaw(void)
{
	MassErrorTol=10;
	MinS_N=2;
	Min_Iso_Num=3;
	Min_Peaks_Num=5;
	GD_CUT=0.01;
	ISO_Win=2.0;
}

IdxRaw::~IdxRaw(void)
{
	
}

void IdxRaw::RemUsed()
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

int IdxRaw::FindIsoCNew(double BaseInt)
{
	size_t pnum=tmpPKL.size();
	if(pnum<=Min_Iso_Num) return NO_PEAKS;
	double oldmz=tmpPKL[0].mz;
	tmpPKL[0].IsUsed=true;
	if(tmpPKL[0].dInt<RINT_CUT*BaseInt)
	{	
		RemUsed();		
		return CPEAK_FAILUR;
	}
	IsoGroup IG;
	size_t i;
	vector<XICPeak> tmpXICPL;
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

	/*if(IG.UList.size()>1)
	{
		printf("\nfind a mutiple charge overlap case:scannum=%d\tbasemz=%lf\n",MSList[CurrentMS].ScanNumber,IG.gPeaks[0].mz);
	}
	*/
	if(IG.gPeaks.size()<Min_Iso_Num)//not a iso group, remove these peaks only
	{
		RemUsed();		
		return CPEAK_FAILUR;
	}
	CurrentBestGD=Fitting(IG,GD_CUT,Min_Iso_Num,tmpXICPL);
	if(CurrentBestGD<GD_CUT) 
	{
		RemUsed();
		return NO_GOOD_CLUSTER;	
	}
	size_t si=tmpXICPL.size();	
	for(i=0;i<si;i++)
	{
		if(tmpXICPL[i].S_N>MinS_N)
		{
			if(BaseInt>DEPSL) tmpXICPL[i].RelInt=tmpXICPL[i].MonoInt/BaseInt;
			AddXICPeaks(tmpXICPL[i]);
		}
	}
	RemUsed();
	return FIND_ONE;
}

int IdxRaw::FindIsoCPmz(double BaseInt,FILE *fp)
{
	size_t pnum=tmpPKL.size();
	if(pnum<Min_Iso_Num) return NO_PEAKS;
	double oldmz=tmpPKL[0].mz;
	tmpPKL[0].IsUsed=true;
	if(tmpPKL[0].dInt<RINT_CUT*BaseInt)
	{	
		RemUsed();		
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

	/*if(IG.UList.size()>1)
	{
		printf("\nfind a mutiple charge overlap case:scannum=%d\tbasemz=%lf\n",MSList[CurrentMS].ScanNumber,IG.gPeaks[0].mz);
	}
	*/
	if(IG.gPeaks.size()<Min_Iso_Num)//not a iso group, remove these peaks only
	{
		RemUsed();		
		return CPEAK_FAILUR;
	}
	CurrentBestGD=Fitting(IG,GD_CUT,Min_Iso_Num,tmpXICPL);
	if(CurrentBestGD<GD_CUT) 
	{
		RemUsed();
		return NO_GOOD_CLUSTER;	
	}
	size_t si=tmpXICPL.size();	
	for(i=0;i<si;i++)
	{
		if(tmpXICPL[i].S_N>MinS_N)
		{
			if(BaseInt>DEPSL) tmpXICPL[i].RelInt=tmpXICPL[i].MonoInt/BaseInt;
			tmpXICPL[i].ScanNum=MSList[CurrentMS].ScanNumber;
			tmpXICPL[i].RT=MSList[CurrentMS].RT;
			tmpXICPL[i].Out2FILE(fp);			
		}
	}
	RemUsed();
	return FIND_ONE;
}

int IdxRaw::FindIsoCPmz()
{
	double BaseInt=0;
	size_t i, pnum=tmpPKL.size();
	if(pnum<Min_Iso_Num) return NO_PEAKS;
	for(i=0;i<pnum;i++)
	{
		if(BaseInt<tmpPKL[i].dInt) BaseInt=tmpPKL[i].dInt;
	}
	double oldmz=tmpPKL[0].mz;
	tmpPKL[0].IsUsed=true;	
	vector<XICPeak> tmpXICPL;
	IsoGroup IG;
	IG.Initial(tmpPKL[0]);
	for(i=1;i<pnum;i++) 
	{
		if(IG.AddPeaks(tmpPKL[i],MassErrorTol))
		{
			tmpPKL[i].IsUsed=true;
			oldmz=tmpPKL[i].mz;
		}	
		else if(tmpPKL[i].mz>oldmz+1.1) break;	
		printf("the unist validate: %d\n",IG.UList.size());
	}
	
	if(IG.gPeaks.size()<Min_Iso_Num)return CPEAK_FAILUR;

	CurrentBestGD=Fitting(IG,GD_CUT,Min_Iso_Num,tmpXICPL);
	if(CurrentBestGD<GD_CUT) return NO_GOOD_CLUSTER;	
	size_t si=tmpXICPL.size();	
	for(i=0;i<si;i++)
	{
		if(BaseInt>DEPSL) tmpXICPL[i].RelInt=tmpXICPL[i].MonoInt/BaseInt;
		tmpXICPL[i].Out2Screen();		
	}	
	return FIND_ONE;
}


int IdxRaw::FindIsoCPmz(int *PreData,int PreNum,FILE *fpout)
{
	double BaseInt=0;
	size_t i, pnum=tmpPKL.size();
	if(pnum<Min_Iso_Num) return NO_PEAKS;
	for(i=0;i<pnum;i++)
	{
		if(BaseInt<tmpPKL[i].dInt) BaseInt=tmpPKL[i].dInt;
	}
	double oldmz=tmpPKL[0].mz;
	tmpPKL[0].IsUsed=true;	
	vector<XICPeak> tmpXICPL;
	IsoGroup IG;
	IG.Initial(tmpPKL[0]);
	for(i=1;i<pnum;i++) 
	{
		if(IG.AddPeaks(tmpPKL[i],MassErrorTol))
		{
			tmpPKL[i].IsUsed=true;
			oldmz=tmpPKL[i].mz;
		}	
		else if(tmpPKL[i].mz>oldmz+1.1) break;	
		printf("the unist validate: %d\n",IG.UList.size());
	}
	
	if(IG.gPeaks.size()<Min_Iso_Num)return CPEAK_FAILUR;

	CurrentBestGD=Fitting(IG,GD_CUT,Min_Iso_Num,tmpXICPL);
	if(CurrentBestGD<GD_CUT) return NO_GOOD_CLUSTER;	
	size_t si=tmpXICPL.size();	
	for(i=0;i<si;i++)
	{
		if(BaseInt>DEPSL) tmpXICPL[i].RelInt=tmpXICPL[i].MonoInt/BaseInt;
		//tmpXICPL[i].Out2Screen();
		//如果存在多个相同电荷的情况，这种方法会得到错误的结果
		//if(PreNum>=tmpXICPL[i].charge) 
		//{			
		//	fprintf(fpout,"%d\t",PreData[tmpXICPL[i].charge-1]);
		//}		
		//else fprintf(fpout,"0\t");
		double tmpDI=1e16;
		int si=-1;
		for(int kk=0;kk<PreNum;kk++)
		{
			double tmpf=fabs(PreData[kk]-tmpXICPL[i].IntU);
			if(tmpf<tmpDI)
			{
				si=kk;
				tmpDI=tmpf;
			}
		}
		if(si!=-1) fprintf(fpout,"%d\t",PreData[si]);
		else fprintf(fpout,"0\t");

		tmpXICPL[i].Out2FileNew(fpout);
	}	
	return FIND_ONE;
}

void IdxRaw::FindMInt(DataPeak *pkl,int pnum)
{
	CurrentBestGD=0;
	tmpPKL.clear();
	double baseInt=GetBasePeak(pkl,pnum);
	Peak pt;	
	for(int i=0;i<pnum;i++)
	{
		pt=pkl[i];
		tmpPKL.push_back(pt);
	}

	int CT=0;
	while(FindIsoCNew(baseInt)!=NO_PEAKS) 
	{
		CT++;
		//printf("\r%d",CT);
	}
	size_t Isomz=_pXICList.size();
	printf("\rScanNum:%ld RT:%lf Peaks:%d XICs:%d Max.GD:%lf",\
		MSList[CurrentMS].ScanNumber,MSList[CurrentMS].RT,CT,Isomz,CurrentBestGD);
}

//not safe enough
int IdxRaw::TFindL(DataPeak *pkl,int b,int e,double fmz)
{
	if(pkl[b].dMass>fmz+0.01) return -1;
	if(pkl[e].dMass<fmz-0.01) return -1;
	if(e-b<=1) return b;	
	int idx=(b+e)/2;
	if(pkl[idx].dMass>=fmz)
	{
		return TFindL(pkl,b,idx,fmz);
	}
	else return TFindL(pkl,idx,e,fmz);
}

int IdxRaw::TFindL(int b,int e,long scan)
{
	if(e-b<=1) return b;	
	int idx=(b+e)/2;
	if(MSList[idx].ScanNumber>=scan)
	{
		return TFindL(b,idx,scan);
	}
	else return TFindL(idx,e,scan);
}


int IdxRaw::FindByPmz(double pmz,DataPeak *pkl,int pnum,FILE *fp)
{
	if(pnum<=0) return 0;
	CurrentBestGD=0;
	tmpPKL.clear();
	double baseInt=GetBasePeak(pkl,pnum);
	int i,targetidx=TFindL(pkl,0,pnum-1,pmz-ISO_Win);
	if(targetidx==-1) return 0;
	Peak pt;
	//extend according the iso pair rules
	for(i=targetidx;i>0;i--)
	{
		double delt=pkl[i].dMass-pkl[i-1].dMass;
		int ch;
		for(ch=MAX_CHARGE;ch>=1;ch--)
		{
			if(fabs(delt*ch-ISO_DIFF[1])<MassErrorTol*pkl[i].dMass/1e6) break;
		}
		if(ch==0) break;
	}
	int begin=i;
	targetidx=TFindL(pkl,0,pnum-1,pmz+ISO_Win);
	if(targetidx==-1) return 0;
	//extend according the iso pair rules
	for(i=targetidx;i<pnum-1;i++)
	{
		double delt=pkl[i+1].dMass-pkl[i].dMass;
		int ch;
		for(ch=MAX_CHARGE;ch>=1;ch--)
		{
			if(fabs(delt*ch-ISO_DIFF[1])<MassErrorTol*pkl[i].dMass/1e6) break;
		}
		if(ch==0) break;
	}
	int end=i;
	for(i=begin;i<=end;i++) 
	{
		pt=pkl[i];
		tmpPKL.push_back(pt);
	}
	int CT=0;	
	int bT=FindIsoCPmz(baseInt,fp);	
	if(bT==FIND_ONE)CT++;
	while(bT!=NO_PEAKS) 
	{	
		bT=FindIsoCPmz(baseInt,fp);	
		if(bT==FIND_ONE)CT++;
		//printf("\r%d",CT);
	}	
	return CT;
}

void IdxRaw::AddXICPeaks(XICPeak &xt,int CH)
{	
	size_t xNum=_pXICList.size();
	size_t k;
	xt.RT=MSList[CurrentMS].RT;
	xt.ScanNum=MSList[CurrentMS].ScanNumber;
	for(k=0;k<xNum;k++)
	{
		if(_pXICList[k].AddPeak(xt,MassErrorTol,CH)) break;
	}
	if(k==xNum)
	{
		XICs xct;
		xct.AddPeak(xt,MassErrorTol,CH);
		_pXICList.push_back(xct);
	}		
}

void IdxRaw::AddXICPeaks(XICPeak &xt)
{	
	size_t xNum=_pXICList.size();
	size_t k;
	xt.RT=MSList[CurrentMS].RT;
	xt.ScanNum=MSList[CurrentMS].ScanNumber;
	for(k=0;k<xNum;k++)
	{
		if(_pXICList[k].AddPeak(xt,MassErrorTol)) break;
	}
	if(k==xNum)
	{
		XICs xct;
		xct.AddPeak(xt,MassErrorTol);
		_pXICList.push_back(xct);
	}		
}

int IdxRaw::FindPP(vector<XICs> *XICList)
{		
	_pXICList.clear();
	size_t i,count=MSList.size();
	if(count<=0) return 0;	
	long PeakNum=0;	
	DataPeak* pDataPeaks = NULL;	
	for(i=0;i<count;i++)
	{		
		PeakNum=MSList[i].PNum;		
		if(PeakNum<=0) continue;
		CurrentMS=i;	
		pDataPeaks=MSList[i].pkl;	
		//if(MSList[CurrentMS].ScanNumber==2443)//debug
		//{//debug
			FindMInt(pDataPeaks,PeakNum);	
		//}//debug
	}

	//CStdioFile debugFile;
	//debugFile.Open("L:\\labelingdata\\BSAT\\debugdata.txt",CFile::modeWrite|CFile::modeCreate|CFile::typeText);

	XICList->clear();
	count=_pXICList.size();
	int CT=0;
	//string tmpStr;
	for(i=0;i<count;i++)
	{
		//_pXICList[i].Output(tmpStr);
		//debugFile.WriteString(tmpStr.c_str());
		CT+=FilterPP(_pXICList[i],XICList);
	}	
	return CT;
}

int IdxRaw::IDXFromScanNum(long ScanNum)
{	
	if(ScanNum<1) return -1;
	size_t MSNum=MSList.size();
	if(MSNum<=0) return -1;
	int e=MSNum-1;
	if(ScanNum>MSList[e].ScanNumber) return -1;
	e=TFindL(0,e,ScanNum);
	if(MSList[e].ScanNumber+10<ScanNum) return -1;//generally, no 10 MS/MS after one MS
	return e;
}

int IdxRaw::ISObyPmz(long ScanNum,double pmz,int ch,FILE *fp)
{		
	_pXICList.clear();
	int idx=IDXFromScanNum(ScanNum);
	if(idx==-1) return 0;

	long PeakNum=0;	
	DataPeak* pDataPeaks = NULL;
	
	PeakNum=MSList[idx].PNum;
	if(PeakNum<=0) return 0;
	pDataPeaks=MSList[idx].pkl;
	CurrentMS=idx;
	//debug for failure
	//if(ScanNum==5857)//5355||ScanNum==5356)//1618,802
	//{
	//	ScanNum=ScanNum;
	//}
	//debug
	fprintf(fp,">The intput is mz=%lf,Scannum=%d, charge=%d\n",pmz,ScanNum,ch);
	int IsFind=FindByPmz(pmz,pDataPeaks,PeakNum,fp);

	printf("\rScanNum:%ld RT:%lf Is Find:%d Max.GD:%lf",\
		MSList[idx].ScanNumber,MSList[idx].RT,IsFind,CurrentBestGD);
	return IsFind;
}

bool IdxRaw::ISObySingleMS(char *MSFile)
{		
	_pXICList.clear();
	FILE *fp;
	fp=fopen(MSFile,"r");
	if(fp==NULL) return false;
	long PeakNum=0;	
	Peak pt;
	pt.base=1;
	pt.noise=1;
	while(!feof(fp))
	{
		int RD=fscanf(fp,"%lf\t%lf\n",&(pt.mz),&(pt.dInt));
		if(RD!=2) continue;
		tmpPKL.push_back(pt);
		
	}
	fclose(fp);
	FindIsoCPmz();	
	return true;
}

bool IdxRaw::ISObySingleMSNew(char *MSFile,char *outfile)
{		
	_pXICList.clear();
	FILE *fp,*fpout;
	fp=fopen(MSFile,"r");
	if(fp==NULL) return false;	
	fpout=fopen(outfile,"w");
	if(fpout==NULL)
	{
		fclose(fp);
		return false;
	}
	long PeakNum=0;	
	Peak pt;
	pt.base=1;
	pt.noise=1;
	char buf[256];
	fgets(buf,255,fp);
	if(strstr(buf,">")!=buf)
	{
		fclose(fp);
		return false;
	}
	int A[3];
	int RD=sscanf(buf,">%d.%d.%d",A,A+1,A+2);
	if(RD!=3)
	{
		fclose(fp);
		return false;
	}
	while(!feof(fp))
	{
		fgets(buf,255,fp);
		if(strstr(buf,">")==buf)
		{
			FindIsoCPmz(A,3,fpout);	
			sscanf(buf,">%d.%d.%d",A,A+1,A+2);	
			tmpPKL.clear();
		}
		else
		{
			RD=sscanf(buf,"%lf\t%lf\n",&(pt.mz),&(pt.dInt));
			if(RD!=2) continue;
			tmpPKL.push_back(pt);
		}		
	}
	fclose(fp);
	fclose(fpout);	
	return true;
}

int IdxRaw::FilterPP(XICs &xct,vector<XICs> *XICList)
{
	vector<XICPeak> XICPL;	
	xct.GetPKL(&XICPL);
	size_t i,count=XICPL.size();
	if(count<=0) return 0;
	XICs splitXIC;
	double RTOld=XICPL[0].RT;
	int CT=0;
	splitXIC.charge=xct.charge;
	splitXIC.AddPeak(XICPL[0]);
	//printf("The current m/z is %lf, charge is %d\n",XICPL[0].mz,xct.charge);
	for(i=1;i<count;i++)
	{		
		double dt=XICPL[i].RT-RTOld;
		if(dt>DT_CUT) 
		{		
			size_t tmpNum=splitXIC.GetPeakNum();	
			//printf("catch a split, the peak number is :%d\n",tmpNum);
			if(tmpNum>Min_Peaks_Num) 
			{			
				XICList->push_back(splitXIC);
				CT++;
			}
			splitXIC.clear();
		}
		splitXIC.AddPeak(XICPL[i]);	
		RTOld=XICPL[i].RT;
	}
	size_t tmpNum=splitXIC.GetPeakNum();	
	//printf("catch a split, the peak number is :%d\n",tmpNum);
	if(tmpNum>Min_Peaks_Num) 
	{
		XICList->push_back(splitXIC);
		CT++;
	}	
	return CT;
}

void IdxRaw::Close()
{
	MSList.clear();
}

double IdxRaw::GetBasePeak(DataPeak *pkl,int m_nArraySize)
{
	double BasePeak=0;
	for(int i=0;i<m_nArraySize;i++)
	{
		if(BasePeak<pkl[i].dIntensity) BasePeak=pkl[i].dIntensity;
	}
	return BasePeak;
}
