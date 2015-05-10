// FTRawFile.cpp: implementation of the FTRawFile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FTRawFile.h"
#include "math.h"
#include "GaussFit.h"
#include "comutil.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FTRawFile::FTRawFile()
{
	charge=0;
	bMass=0;
	sFName[0]='\0';
	MET=10;//default  is 10PPM
	fp_cal_data=NULL;
	IsOutput=false;

}

FTRawFile::~FTRawFile()
{
	if(IsOutput)
	{
		if(fp_cal_data!=NULL) fclose(fp_cal_data);		
	}
}

bool FTRawFile::Initial(CString sRwafile)
{
	if(Open(sRwafile)==0)
	{
		long lRet = -1;
		GetNumberOfControllers(&lRet);
		if(lRet)
		{
			SetCurrentController(0,1);
		}
		else
		{
			Close();
			return false;
		}
	}
	else return false;
	SNFirst=0;
	GetFirstSpectrumNumber(&SNFirst);	
	SNLast=0;	
	GetLastSpectrumNumber(&SNLast);

	if(SNLast<=SNFirst) 
	{
		Close();
		return false;
	}
	strcpy(sFName,sRwafile.GetBuffer());
	sRwafile.ReleaseBuffer();
	//remove the extensition
	sFName[sRwafile.GetLength()-4]='\0';
	return true;
}

double FTRawFile::match(mPeak ISOP[MAX_ISO],double IsoDisT[MAX_ISO])
{
	double sumx,sumy,sumxy,sumI;
	sumx=0;
	sumy=0;
	sumxy=0;
	sumI=0;
	int i;	
	for(i=0;i<MAX_ISO;i++)
	{
		sumI+=ISOP[i].dIntensity;
	}
	if(sumI<1e-6) return 0;
	for(i=0;i<MAX_ISO;i++)
	{
		//if(ISOP[i].dIntensity<=0) break;
		double NormalI=ISOP[i].dIntensity/sumI;
		sumx+=NormalI*NormalI;
		sumy+=IsoDisT[i]*IsoDisT[i];
		sumxy+=NormalI*IsoDisT[i];
	}	
	//if(i<ISO_CUT) return 0;
	sumx*=sumy;
	sumx=sqrt(sumx);
	if(sumx<=0) return 0;	
	return sumxy/sumx;
}

double FTRawFile::PreFindMInt(IsoCluster &it)
{
	if(charge<=0) return 0;	
	double dm=0.05;//charge*bMass*MET/1e6;//20ppm
	double step=1.007825/charge;
	double subLP[MAX_ISO+3],subUP[MAX_ISO+3];
	mPeak tmpISO[MAX_ISO+3];
	long i,k;
	for(i=0;i<MAX_ISO+3;i++)
	{
		double Cen=bMass+(i-3)*step;
		subLP[i]=Cen-dm;
		subUP[i]=Cen+dm;
	}
	
	//double basePeak=0;
	size_t size=MSSpec.size();
	//for(i=0;i<size;i++) 
	//{
	//	if(basePeak<MSSpec[i].dIntensity) basePeak=MSSpec[i].dIntensity;
	//}
	int Begin=TFindL(MSSpec,0,size,subLP[0]);
	int End=TFindL(MSSpec,Begin,size,subUP[MAX_ISO-1]);
	if(End+1<size) End++;
	for(i=Begin;i<End;i++)
	{		
		for(k=0;k<MAX_ISO+3;k++)
		{
			if(MSSpec[i].dMass>subLP[k]&&MSSpec[i].dMass<subUP[k])
			{
				if(tmpISO[k].dIntensity<MSSpec[i].dIntensity) 
				{
					tmpISO[k]=MSSpec[i];	
				}
			}
		}			
	}

	int basei=3;
	//double IsoDisT[MAX_ISO];
	CalDefIsoDis(bMass*charge,IsoDisT);
	double gd_max=match(tmpISO+3,IsoDisT);;
	double gd;
	//find the mono peaks
	for(i=2;i>=0;i--)
	{
		gd=match(tmpISO+i,IsoDisT);
		if(gd>gd_max) 
		{
			basei=i;
			gd_max=gd;
		}
	}
	if(gd_max<GD_CUT) return 0;
	for(i=0;i<MAX_ISO;i++)
	{
		it.ISOP[i]=tmpISO[i+basei];
	}

/*	if(basePeak>1e-4)
	{
		for(i=0;i<MAX_ISO;i++)
		{
			it.ISOP[i].rInt=it.ISOP[i].dIntensity/basePeak;
		}
	}*/	
	if(it.ISOP[0].rInt<RINT_CUT) return 0;
	it.charge=charge;
	bMass=it.ISOP[0].dMass;//record the refined bMass
	if(bMass<=MASS_LIM) return 0;
	return gd_max;
}

double FTRawFile::PreFindMIntSim()
{	
	double dm=0.05;//the float point data reduction error range	
	size_t size=MSSpec.size();
	int Begin=TFindL(MSSpec,0,size,bMass-dm);
	int End=TFindL(MSSpec,Begin,size,bMass+dm);
	if(End+1<size) End++;
	double tmpMaxInt=0;
	int tmpIDX=-1;
	for(int i=Begin;i<End;i++)
	{
		if(tmpMaxInt<MSSpec[i].dIntensity) 
		{
			tmpMaxInt=MSSpec[i].dIntensity;
			tmpIDX=i;
		}
	}
	if(tmpIDX==-1) return 0;
	return MSSpec[tmpIDX].dMass;
}

bool FTRawFile::FindMInt(IsoCluster &it)
{
	if(charge<=0) return false;	
	double dm=charge*bMass*MET/1e6;//20ppm
	double step=1.007825/charge;
	double subLP[MAX_ISO],subUP[MAX_ISO];
	long i,k;
	for(i=0;i<MAX_ISO;i++)
	{
		it.ISOP[i].dIntensity=0;
		it.ISOP[i].rInt=0;
		it.ISOP[i].dMass=0;
		double Cen=bMass+i*step;
		subLP[i]=Cen-dm;
		subUP[i]=Cen+dm;
	}	
	//double basePeak=0;
	long size=MSSpec.size();
	//for(i=0;i<size;i++) 
	//{
	//	if(basePeak<MSSpec[i].dIntensity) basePeak=MSSpec[i].dIntensity;
	//}

	int Begin=TFindL(MSSpec,0,size,subLP[0]);
	int End=TFindL(MSSpec,Begin,size,subUP[MAX_ISO-1]);
	if(End+1<size) End++;
	for(i=Begin;i<End;i++)
	{		
		for(k=0;k<MAX_ISO;k++)
		{
			if(MSSpec[i].dMass>subLP[k]&&MSSpec[i].dMass<subUP[k])
			{
				if(it.ISOP[k].dIntensity<MSSpec[i].dIntensity) 
				{
					it.ISOP[k]=MSSpec[i];	
				}
			}
		}			
	}
	//if(basePeak>1e-4)
	//{
	//	for(i=0;i<MAX_ISO;i++)
	//	{
	//		it.ISOP[i].rInt=it.ISOP[i].dIntensity/basePeak;
	//	}
	//}
	if(it.ISOP[0].rInt<RINT_CUT) return false;
	double IsoDisT[MAX_ISO];
	CalDefIsoDis(bMass*charge,IsoDisT);
	if(it.Match(IsoDisT)<GD_CUT) return false;
	it.charge=charge;
	return true;
}

bool FTRawFile::FindMIntFit(IsoCluster &it)
{
	if(charge<=0) return false;	
	double dm=charge*bMass*MET/1e6;//20ppm
	double step=1.007825/charge;
	double subLP[MAX_ISO],subUP[MAX_ISO];
	long i,k;
	for(i=0;i<MAX_ISO;i++)
	{
		it.ISOP[i].dIntensity=0;
		it.ISOP[i].rInt=0;
		it.ISOP[i].dMass=0;
		double Cen=bMass+i*step;
		subLP[i]=Cen-dm;
		subUP[i]=Cen+dm;
	}	
	//double basePeak=0;
	long size=MSSpec.size();
	//for(i=0;i<size;i++) 
	//{
	//	if(basePeak<MSSpec[i].dIntensity) basePeak=MSSpec[i].dIntensity;
	//}

	int Begin=TFindL(MSSpec,0,size,subLP[0]);
	int End=TFindL(MSSpec,Begin,size,subUP[MAX_ISO-1]);
	if(End+1<size) End++;
	for(i=Begin;i<End;i++)
	{		
		for(k=0;k<MAX_ISO;k++)
		{
			if(MSSpec[i].dMass>subLP[k]&&MSSpec[i].dMass<subUP[k])
			{
				if(it.ISOP[k].dIntensity<MSSpec[i].dIntensity) 
				{
					it.ISOP[k]=MSSpec[i];	
				}
			}
		}			
	}
	//if(basePeak>1e-4)
	//{
	//	for(i=0;i<MAX_ISO;i++)
	//	{
	//		it.ISOP[i].rInt=it.ISOP[i].dIntensity/basePeak;
	//	}
	//}
	if(it.ISOP[0].rInt<RINT_CUT) return false;
	double IsoDisT[MAX_ISO];
	CalDefIsoDis(bMass*charge,IsoDisT);
	if(it.Match(IsoDisT)<GD_CUT) return false;
	it.charge=charge;
	return true;
}

int FTRawFile::TFindL(vector<mPeak> &pkl,int b,int e,double fmz)
{
	if(e-b<=1) return b;	
	int idx=(b+e)/2;
	if(pkl[idx].dMass>=fmz)
	{
		return TFindL(pkl,b,idx,fmz);
	}
	else return TFindL(pkl,idx,e,fmz);
}

bool FTRawFile::RefinePISO(IsoCluster &it)
{
	if(FindMInt(it)) return true;
	return false;
}

bool FTRawFile::PreRefinePISO(IsoCluster &it)
{
	double gd_tmp,gd_max=0.5;
	int ch_final=0;
	IsoCluster it_tmp;
	double tmpBmass=bMass;
	for(charge=MAX_CH;charge>=1;charge--)
	{			
		gd_tmp=PreFindMInt(it);
		if(gd_tmp>gd_max) 
		{
			ch_final=charge;
			gd_max=gd_tmp;
			it_tmp=it;
		}
		bMass=tmpBmass;//reset the bMass
	}
	if(ch_final==0) return false;
	charge=ch_final;
	it=it_tmp;
	bMass=it.ISOP[0].dMass;
	return true;
}

void FTRawFile::GetIsoCData(IsoCluster *it,vector<double> *dMass,double C[5])
{
	//MAX_ISO,only use the Mono mass	
	for(int i=0;i<MAX_ISO;i++)
	{
		if(it->ISOP[i].rInt<1e-6) break;
		double mze=it->ISOP[i].dMass/MZE_SCALE;
		double tic=it->pTIC/TIC_SCALE;
		double rt=it->RT;
		double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*rt-ISO_DIFF[i]/(it->charge);
		dMass->push_back(cmass);
		if(IsOutput)
		{
			fprintf(fp_cal_data,"%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",it->ScanNum,it->charge,it->RT,it->pTIC,it->ISOP[i].rInt,it->ISOP[i].dIntensity,it->ISOP[i].dMass,cmass-bMass,it->ISOP[i].dMass-bMass);
		}
	}
}

void FTRawFile::GetIsoCDataExt(IsoCluster *it,vector<double> *dMass,double C[pdNum])
{
	//MAX_ISO,only use the Mono mass	
	int i=0;
	double commoncmass=0;	
	double tic=it->pTIC/TIC_SCALE;
	double rt=sqrt(it->RT);
	for(i=0;i<MAX_PAR_NUM;i++) commoncmass+=C[i+5]*it->FTStatus[i]+C[LIN_PAR_NUM+1]*it->GetIsoNum()+C[4]*rt;
	for(i=0;i<MAX_ISO;i++)
	{
		if(it->ISOP[i].rInt<1e-2) break;
		double mze=it->ISOP[i].dMass/MZE_SCALE;	
		double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic-ISO_DIFF[i]/(it->charge)+ \
			commoncmass+C[LIN_PAR_NUM]*it->ISOP[i].rInt;		
		dMass->push_back(cmass);
		if(IsOutput)
		{
			fprintf(fp_cal_data,"%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",it->ScanNum,it->charge,it->RT,it->pTIC,it->ISOP[i].rInt,it->ISOP[i].dIntensity,it->ISOP[i].dMass,cmass-bMass,it->ISOP[i].dMass-bMass);
		}
	}
}

void FTRawFile::GetXICData(vector<IsoCluster> *XICList,vector<double> *dMass,double C[5])
{
	dMass->clear();
	size_t i,s=XICList->size();	
	for(i=0;i<s;i++)
	{
		GetIsoCData(&(XICList->at(i)),dMass,C);
	}
}

void FTRawFile::GetXICDataExt(vector<IsoCluster> *XICList,vector<double> *dMass,double C[pdNum])
{
	dMass->clear();
	size_t i,s=XICList->size();	
	for(i=0;i<s;i++)
	{
		GetIsoCDataExt(&(XICList->at(i)),dMass,C);
	}
}


void FTRawFile::CalMST(vector<double> *dMass,double *m,double *s)
{
	*m=0;
	*s=0;
	size_t i,size=dMass->size();
	if(size<=0) return;
	for(i=0;i<size;i++)
	{
		*m+=dMass->at(i)/1000;
	}
	*m/=size;

	for(i=0;i<size;i++)
	{
		double t=dMass->at(i)/1000-*m;
		*s+=t*t;
	}
	*s=sqrt(*s/size);
	*m=*m*1000;
	*s=*s*1000;
}

double FTRawFile::w_mean(vector<double> *x, char *W)
{
   int num=0;
   double m=0;
   size_t n=x->size();
   for(size_t i=0;i<n;i++)
   {
	   if(W[i])
	   {
		   m+=x->at(i);
		   num++;
	   }
   }
   if(num>0) m/=num;
   return m;
}

double FTRawFile::w_std(vector<double> *x,char *W,double m)
{
	double std=0;
	int num=0;
	size_t n=x->size();
	for(size_t i=0;i<n;i++)
	{
	   if(W[i])
	   {
		   double t=x->at(i)-m;
		   std+=t*t;
		   num++;
	   }
	}
	if(num>0) std=sqrt(std/num);
	return std;
}

void FTRawFile::UpdateResW(vector<double> *x,char *W,double mean,double sigma)
{
	size_t n=x->size();
	for(size_t i=0;i<n;i++)
	{
		if(fabs(x->at(i)-mean)>STD_TUNE*sigma) W[i]=0;
		else W[i]=1;
	}
}

void FTRawFile::rCalMST(vector<double> *dMass,double *m,double *s)
{
	size_t n=dMass->size();
	if(n<=0) 
	{
		*m=0;
		*s=0;
		return;
	}
	char *w;
	w=new char[n];
	size_t i;
	for(i=0;i<n;i++) w[i]=1;
	double mean=w_mean(dMass,w);
	double sigma=w_std(dMass,w,mean);
	double tol=1;
	double sigmaold=sigma;
	int count=0;
	while(tol>1e-3)
	{	
		UpdateResW(dMass,w,mean,sigma);
		mean=w_mean(dMass,w);
		sigma=w_std(dMass,w,mean);
		if(sigmaold<=0) sigmaold=STD_EPS;
		tol=fabs(sigma-sigmaold)/sigmaold;
		sigmaold=sigma;
		count++;
		if(count>100) break;
	}
	*m=mean;
	*s=sigma;	
	delete []w;
}

double FTRawFile::GetMSTIC()
{
	double pTIC=0;
	size_t size=MSSpec.size();
	for(int i=0;i<size;i++) pTIC+=MSSpec[i].dIntensity;
	return pTIC;
}

int FTRawFile::MSLevel(long ScanNum)
{
	CString m_sActualFilter;
	BSTR bstrFilter = NULL;
	GetFilterForScanNum(ScanNum, &bstrFilter);
	m_sActualFilter = bstrFilter;
	SysFreeString(bstrFilter);
	bstrFilter=NULL;
	if(m_sActualFilter.Find("Full ms [")!=-1) return 1;
	if(m_sActualFilter.Find("Full ms2 ")!=-1) return 2;		
	if(m_sActualFilter.Find("Full ms3 ")!=-1) return 3;	
	return 0;
}

bool FTRawFile::GetMS2pMass(long ScanNum)
{
	CString m_sActualFilter;
	BSTR bstrFilter = NULL;
	GetFilterForScanNum(ScanNum, &bstrFilter);
	m_sActualFilter = bstrFilter;
	SysFreeString(bstrFilter);
	bstrFilter=NULL;
	int idx=m_sActualFilter.Find("Full ms2 ");
	if(idx==-1)  return false;
	idx+=9;
	m_sActualFilter=m_sActualFilter.Mid(idx);
	if(sscanf(m_sActualFilter.GetBuffer(),"%lf",&bMass)!=1)
	{
		m_sActualFilter.ReleaseBuffer();
		return false;
	}
	m_sActualFilter.ReleaseBuffer();
	return true;
}

bool FTRawFile::ReadMSLabel(long ScanNum,vector<mPeak> &MSSpecTmp)
{	
	if(ScanNum>SNLast||ScanNum<SNFirst) return false;
	MSSpecTmp.clear();
	mPeak pt;
	long	nScanNumber = ScanNum; // get  the label data of the first scan.
	int		dim, inx, charge;
	double		*pdval;
	unsigned char	*pcval;
	SAFEARRAY   	*parray, *parray2;
	_variant_t	vSpecData, vFlags;
	VARIANT	varLabels, *pvarLabels;
	VARIANT	varFlags, *pvarFlags;
		
	double		dMass, dInt;
	unsigned char	cMerged, cFragmented, cReference, cException, cModified, cSaturated;
	TCHAR		flags[7];
	float		fRes, fBase, fNoise;
			
	pvarLabels = &varLabels;
	pvarFlags	= &varFlags;
	GetLabelData(pvarLabels, pvarFlags, &nScanNumber);

	vSpecData	= pvarLabels;
	parray		= vSpecData.parray;
	dim		= parray->rgsabound[0].cElements;	
	pdval		= (double *) parray->pvData;
			
	if(pvarFlags)
	{
		vFlags	= pvarFlags;
		parray2	= vFlags.parray;
		pcval	= (unsigned char *) parray2->pvData;
	}		
			
	for (inx = 0; inx < dim; inx++)
	{
		dMass		=  (double)	pdval[((inx)*6)+0] ;
		dInt		=  (double)	pdval[((inx)*6)+1] ;
		fRes		=  (float) 	pdval[((inx)*6)+2] ;
		fBase		=  (float) 	pdval[((inx)*6)+3] ;
		fNoise		=  (float) 	pdval[((inx)*6)+4] ;
		charge		=  (int)		pdval[((inx)*6)+5] ;
		if(pvarFlags)
		{
			cSaturated	=  (unsigned char)  pcval[((inx)*6)+0] ;
			cFragmented	=  (unsigned char)  pcval[((inx)*6)+1] ;
			cMerged	=  (unsigned char)  pcval[((inx)*6)+2] ;
			cException	=  (unsigned char)  pcval[((inx)*6)+3] ;
			cReference	=  (unsigned char)  pcval[((inx)*6)+4] ;
			cModified	=  (unsigned char)  pcval[((inx)*6)+5] ;
						
			// write the flags into a String
			flags[0] = _T('\0');
			if(cSaturated)	_tcscat(flags, _T("S"));
			if(cFragmented)	_tcscat(flags, _T("F"));
			if(cMerged)	_tcscat(flags, _T("M"));
			if(cException)_tcscat(flags, _T("E"));
			if(cReference)_tcscat(flags, _T("R"));
			if(cModified)_tcscat(flags, _T("O"));
			//printf("%lf\t%lf\t%f\t%f\t%f\t%d\t%s\n",dMass,dInt,fRes,fBase,fNoise,charge,flags);
			pt.dMass=dMass;
			pt.dIntensity=dInt;
			pt.baseLine=fBase;	
			MSSpecTmp.push_back(pt);
		}
	}
	return true;
}

bool FTRawFile::ReadMS(long ScanNum,vector<mPeak> &MSSpecTmp)
{
	if(ScanNum>SNLast||ScanNum<SNFirst) return false;
	VARIANT m_varMassList;		
	VARIANT varPeakFlags;
	int m_nThresholdType=0;
    UINT m_nThresholdValue=0;
	long m_nArraySize=0;
    long m_nMaxNumPeaks=0;
	bool m_bCentroidScanData=false;	
	double pdCentroidPeakWidth=0;
	VariantInit(&m_varMassList);
	VariantInit(&varPeakFlags);	
	DataPeak* pDataPeaks;
	pDataPeaks = NULL;
	//read the MS scan
	GetMassListFromScanNum(&ScanNum, 
		NULL, 
		m_nThresholdType, 
		m_nThresholdValue, 
		m_nMaxNumPeaks,
		m_bCentroidScanData,
		&pdCentroidPeakWidth,
		&m_varMassList, 
		&varPeakFlags,
		&m_nArraySize);

	if(m_nArraySize<=0) return false;
	SAFEARRAY FAR* psa=m_varMassList.parray;	
	SafeArrayAccessData(psa,(void**)(&pDataPeaks));
	//refine the pmz	
	MSSpecTmp.clear();
	mPeak pt;
	double MaxInt=0;
	int i;
	for(i=0;i<m_nArraySize;i++)
	{
		if(pDataPeaks[i].dIntensity<=0.01) continue;
		if(MaxInt<pDataPeaks[i].dIntensity) MaxInt=pDataPeaks[i].dIntensity;
		pt=pDataPeaks[i];
		MSSpecTmp.push_back(pt);
	}
	m_nArraySize=MSSpecTmp.size();
	for(i=0;i<m_nArraySize;i++) MSSpecTmp[i].rInt=MSSpecTmp[i].dIntensity/MaxInt;
	SafeArrayUnaccessData(psa);
	if(m_varMassList.vt!=VT_EMPTY)
	{
		psa = m_varMassList.parray;
		m_varMassList.parray = NULL;
		// Delete the SafeArray
		SafeArrayDestroy(psa);
	}
	if(varPeakFlags.vt!=VT_EMPTY)
	{
		psa=varPeakFlags.parray;
		varPeakFlags.parray = NULL;
		// Delete the SafeArray
		SafeArrayDestroy(psa);
	}
	return true;
}

double FTRawFile::GetPMassL(long ScanNum,double mzeP)
{
	vector<mPeak> MSSpecTmp;
	if(MSLevel(ScanNum)!=1) return 0;	
	if(!ReadMSLabel(ScanNum,MSSpecTmp)) return 0;
	size_t st=MSSpecTmp.size();
	double dm=0.1;
	double breturn=0;
	for(size_t i=0;i<st;i++)//find the most closeable one
	{
		double dtm=fabs(mzeP-MSSpecTmp[i].dMass);
		if(dtm<dm)
		{
			breturn=MSSpecTmp[i].dMass;
			dm=dtm;
		}
	}
	return breturn;
}

double FTRawFile::GetGCError(double C[5], long ScanNum,double pMZ,int ch)
{
	if(MSLevel(ScanNum)!=2) return 0;	
	//refine the parent mz value in the MS scan
	bMass=pMZ;
	charge=ch;
	CalDefIsoDis(pMZ*ch,IsoDisT);
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return 0;	
	if(!ReadMS(tempSNo,MSSpec)) return 0;	
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;	
	////////////////////////////
	if(!FindMIntFtISO(it)) return 0;
	it.pTIC=GetMSTIC();	
	//it.ScanNum=1;//to record is MS/MS precusor
	it.ScanNum=tempSNo;
	it.charge=ch;

	RTFromScanNum(tempSNo,&(it.RT));
	XICList_F.push_back(it);
	//it.ScanNum=0;//to record is not MS/MS precusor
	int Interupt_times=0;
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);		
		if(FindMIntFtISO(it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			it.ScanNum=tempSNo;
			if(it.goodness>0.9)	XICList_F.push_back(it);
			Interupt_times=0;			
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
    XICList.push_back(XICList_F[0]);
	Interupt_times=0;
	while(tempSNo<SNLast)
	{
		tempSNo++;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);				
		if(FindMIntFtISO(it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			it.ScanNum=tempSNo;
			if(it.goodness>0.9)	XICList.push_back(it);
			Interupt_times=0;	
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	vector<double> dMass;

	GetXICData(&XICList,&dMass,C);
	double m,s;
	rCalMST(&dMass,&m,&s);	
	if(IsOutput)
	{
		fprintf(fp_cal_data,">%d\t%d\t%lf\t%lf\t%lf\n",dMass.size(),ch,pMZ,(m-pMZ)*1e6/pMZ,3*s);
	}
	return (m-pMZ)*1e6/pMZ;
}


double FTRawFile::GetGCErrorExt(double C[pdNum], long ScanNum,double pMZ,int ch)
{
	if(MSLevel(ScanNum)!=2) return 0;	
	//refine the parent mz value in the MS scan
	bMass=pMZ;
	charge=ch;
	CalDefIsoDis(pMZ*ch,IsoDisT);
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return 0;	
	if(!ReadMS(tempSNo,MSSpec)) return 0;	
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;	
	////////////////////////////
	if(!FindMIntFtISO(it)) return 0;
	if(!GetFTingData(tempSNo,it)) return 0;
	//it.pTIC=GetMSTIC();	
	//it.ScanNum=1;//to record is MS/MS precusor
	it.ScanNum=tempSNo;
	it.charge=ch;

	//RTFromScanNum(tempSNo,&(it.RT));
	XICList_F.push_back(it);
	//it.ScanNum=0;//to record is not MS/MS precusor
	int Interupt_times=0;
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);		
		if(FindMIntFtISO(it)&&GetFTingData(tempSNo,it)) 
		{
			//it.pTIC=GetMSTIC();
			//RTFromScanNum(tempSNo,&(it.RT));
			it.ScanNum=tempSNo;
			if(it.goodness>0.9)	XICList_F.push_back(it);
			Interupt_times=0;			
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
    XICList.push_back(XICList_F[0]);
	Interupt_times=0;
	while(tempSNo<SNLast)
	{
		tempSNo++;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);				
		if(FindMIntFtISO(it)&&GetFTingData(tempSNo,it)) 
		{
			//it.pTIC=GetMSTIC();
			//RTFromScanNum(tempSNo,&(it.RT));
			it.ScanNum=tempSNo;
			if(it.goodness>0.9)	XICList.push_back(it);
			Interupt_times=0;	
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	vector<double> dMass;

	//debug  2010.10.27
	//if(XICList.size()==0)
	//{
	//	//to find why?
	//	int pp=0;
	//	pp++;
	//}
	//end debug

	GetXICDataExt(&XICList,&dMass,C);
	double m,s;
	rCalMST(&dMass,&m,&s);	
	if(IsOutput)
	{
		fprintf(fp_cal_data,">%d\t%d\t%lf\t%lf\t%lf\n",dMass.size(),ch,pMZ,(m-pMZ)*1e6/pMZ,3*s);
	}
	return (m-pMZ)*1e6/pMZ;
}


//added on 2010.11.10
bool FTRawFile::GetSkewData(long ScanNum,double pMZ)
{
	if(MSLevel(ScanNum)!=2) return 0;	
	//refine the parent mz value in the MS scan
	bMass=pMZ;	
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;	
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	mPeak mt;
	if(!FindMIntFt(mt)) return false;
	if(IsOutput)
	{
		fprintf(fp_cal_data,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",ScanNum,pMZ,mt.RSquare,mt.rInt,mt.skew[0],mt.skew[1],mt.skew[2],mt.skew[3]);
	}
	return true;
}


bool FTRawFile::SCalibrate(double C[5],FILE *MGFfp,long ScanNum)
{		
	if(!GetMS2pMass(ScanNum)) return false;	
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;
	if(!ReadMS(tempSNo,MSSpec)) return false;
	//refine the pmz
	IsoCluster it;	
	if(!PreRefinePISO(it)) return false;	
	it.pTIC=GetMSTIC();
	RTFromScanNum(tempSNo,&(it.RT));
	double mze=it.ISOP[0].dMass/MZE_SCALE;
	double tic=it.pTIC/TIC_SCALE;
	double rt=it.RT;
	double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*rt;	
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d\n",sFName,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",cmass*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",cmass);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);
	
	if(ReadMS(ScanNum,MSSpec))
	{
		size_t size=MSSpec.size();
		for(int i=0;i<size;i++)
		{
			if(MSSpec[i].dIntensity<1e-2) continue;
			fprintf(MGFfp,"%lf\t%lf\n",MSSpec[i].dMass,MSSpec[i].dIntensity);
		}

	}	
	fprintf(MGFfp,"END IONS\n");
	return true;
}


bool FTRawFile::SCalibrateExt(double C[pdNum],FILE *MGFfp,long ScanNum)
{		
	if(!GetMS2pMass(ScanNum)) return false;	
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;
	if(!ReadMS(tempSNo,MSSpec)) return false;
	//refine the pmz
	IsoCluster it;	
	if(!PreRefinePISO(it)) return false;
	if(!GetFTingData(tempSNo,it)) return false;

	//it.pTIC=GetMSTIC();
	//RTFromScanNum(tempSNo,&(it.RT));
	double mze=it.ISOP[0].dMass/MZE_SCALE;
	double tic=it.pTIC/TIC_SCALE;
	//double rt=it.RT;
	double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*it.RT+C[LIN_PAR_NUM]*it.ISOP[0].rInt+C[LIN_PAR_NUM+1]*it.GetIsoNum();
	for(int i=0;i<MAX_PAR_NUM;i++) cmass+=C[i+5]*it.FTStatus[i];
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d\n",sFName,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",cmass*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",cmass);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);
	
	if(ReadMS(ScanNum,MSSpec))
	{
		size_t size=MSSpec.size();
		for(int i=0;i<size;i++)
		{
			if(MSSpec[i].dIntensity<1e-2) continue;
			fprintf(MGFfp,"%lf\t%lf\n",MSSpec[i].dMass,MSSpec[i].dIntensity);
		}

	}	
	fprintf(MGFfp,"END IONS\n");
	return true;
}

double FTRawFile::SCalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI)
{		
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return ERROR_MASS;
	if(!ReadMS(tempSNo,MSSpec)) return ERROR_MASS;
	//refine the pmz
	bMass=mzI;
	charge=ch;
	IsoCluster it;	
	if(PreFindMInt(it)<0.8) return ERROR_MASS;
	if(!GetFTingData(tempSNo,it)) return ERROR_MASS;

	//it.pTIC=GetMSTIC();
	//RTFromScanNum(tempSNo,&(it.RT));
	double mze=it.ISOP[0].dMass/MZE_SCALE;
	double tic=it.pTIC/TIC_SCALE;
	//double rt=it.RT;
	double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*it.RT+C[LIN_PAR_NUM]*it.ISOP[0].rInt+C[LIN_PAR_NUM+1]*it.GetIsoNum();
	for(int i=0;i<MAX_PAR_NUM;i++) cmass+=C[i+5]*it.FTStatus[i];	
	return cmass;
}

bool FTRawFile::ECalibrate(double C[5],FILE *MGFfp,long ScanNum)
{
	if(!GetMS2pMass(ScanNum)) return false;	
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;
	if(!ReadMS(tempSNo,MSSpec)) return false;
	//refine the pmz	
	double bmass=PreFindMIntSim();
	double pTIC=GetMSTIC();
	double rt;
	RTFromScanNum(tempSNo,&(rt));
	double mze=bmass/MZE_SCALE;
	double tic=pTIC/TIC_SCALE;	
	double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*rt;	
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d\n",sFName,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",cmass*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",cmass);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);
	
	if(ReadMS(ScanNum,MSSpec))
	{
		size_t size=MSSpec.size();
		for(int i=0;i<size;i++)
		{
			if(MSSpec[i].dIntensity<1e-2) continue;
			fprintf(MGFfp,"%lf\t%lf\n",MSSpec[i].dMass,MSSpec[i].dIntensity);
		}

	}	
	fprintf(MGFfp,"END IONS\n");
	return true;
}

bool FTRawFile::ECalibrateExt(double C[pdNum],FILE *MGFfp,long ScanNum)
{
	if(!GetMS2pMass(ScanNum)) return false;	
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;
	if(!ReadMS(tempSNo,MSSpec)) return false;
	//refine the pmz	
	double bmass=PreFindMIntSim();
	IsoCluster it;	
	if(!GetFTingData(tempSNo,it)) return false;
	//double pTIC=GetMSTIC();
	//double rt;
	//RTFromScanNum(tempSNo,&(rt));
	double mze=bmass/MZE_SCALE;
	double tic=it.pTIC/TIC_SCALE;	
	double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*it.RT+C[LIN_PAR_NUM]*it.ISOP[0].rInt+C[LIN_PAR_NUM+1]*it.GetIsoNum();	
	for(int i=0;i<MAX_PAR_NUM;i++) cmass+=C[i+5]*it.FTStatus[i];
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d\n",sFName,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",cmass*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",cmass);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);
	
	if(ReadMS(ScanNum,MSSpec))
	{
		size_t size=MSSpec.size();
		for(int i=0;i<size;i++)
		{
			if(MSSpec[i].dIntensity<1e-2) continue;
			fprintf(MGFfp,"%lf\t%lf\n",MSSpec[i].dMass,MSSpec[i].dIntensity);
		}

	}	
	fprintf(MGFfp,"END IONS\n");
	return true;
}

double FTRawFile::ECalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI)
{	
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return ERROR_MASS;
	if(!ReadMS(tempSNo,MSSpec)) return ERROR_MASS;
	//refine the pmz	
	bMass=mzI;
	charge=ch;
	IsoCluster it;	
	if(!GetFTingData(tempSNo,it)) return ERROR_MASS;
	double pTIC=GetMSTIC();
	double rt;
	RTFromScanNum(tempSNo,&(rt));
	double mze=mzI/MZE_SCALE;
	double tic=pTIC/TIC_SCALE;	
	double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*rt+C[LIN_PAR_NUM]*it.ISOP[0].rInt+C[LIN_PAR_NUM+1]*it.GetIsoNum();;	
	for(int i=0;i<MAX_PAR_NUM;i++) cmass+=C[i+5]*it.FTStatus[i];
	return cmass;	
}


bool FTRawFile::XICCalibrateBack(double C[5],FILE *MGFfp,long ScanNum)
{
	if(!GetMS2pMass(ScanNum)) return false;		
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;
	while(!PreRefinePISO(it)) 
	{
		tempSNo--;
		if(tempSNo<SNFirst) return false;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);	
	}
	it.pTIC=GetMSTIC();	
	RTFromScanNum(tempSNo,&(it.RT));
	XICList_F.push_back(it);
	int Interupt_times=0;
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);			
		if(RefinePISO(it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			XICList_F.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
	XICList.push_back(XICList_F[0]);
	Interupt_times=0;
	while(tempSNo<SNLast)
	{
		tempSNo++;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);		
		if(RefinePISO(it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			XICList.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;	
		if(Interupt_times>=IT_LIM) break;
	}

	if(XICList.size()<=0) return false;

	vector<double> dMass;
	GetXICData(&XICList,&dMass,C);
	double m,s;
	CalMST(&dMass,&m,&s);
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d_%lf.%d.%d\n",sFName,dMass.size(),s,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",m*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",m);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);
	
	if(ReadMS(ScanNum,MSSpec))
	{
		size=MSSpec.size();
		for(int i=0;i<size;i++)
		{
			if(MSSpec[i].dIntensity<1e-2) continue;
			fprintf(MGFfp,"%lf\t%lf\n",MSSpec[i].dMass,MSSpec[i].dIntensity);
		}
	}		
	fprintf(MGFfp,"END IONS\n");
	return true;
}

bool FTRawFile::XICCalibrate(double C[5],FILE *MGFfp,long ScanNum)
{
	if(!GetMS2pMass(ScanNum)) return false;		
    //refine the parent mz value in the MS scan	
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;
	while(!PreFindMIntFtISO(it)) 
	{
		tempSNo--;
		if(tempSNo<SNFirst) return false;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);	
	}
	it.pTIC=GetMSTIC();	
	RTFromScanNum(tempSNo,&(it.RT));
	XICList_F.push_back(it);
	int Interupt_times=0;
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);			
		if(FindMIntFtISO(it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			XICList_F.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
	XICList.push_back(XICList_F[0]);
	Interupt_times=0;
	while(tempSNo<SNLast)
	{
		tempSNo++;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);		
		if(FindMIntFtISO(it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			XICList.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;	
		if(Interupt_times>=IT_LIM) break;
	}

	if(XICList.size()<=0) return false;

	vector<double> dMass;
	GetXICData(&XICList,&dMass,C);
	double m,s;
	CalMST(&dMass,&m,&s);
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d_%lf.%d.%d\n",sFName,dMass.size(),s,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",m*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",m);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);
	
	if(ReadMS(ScanNum,MSSpec))
	{
		size=MSSpec.size();
		for(int i=0;i<size;i++)
		{
			if(MSSpec[i].dIntensity<1e-2) continue;
			fprintf(MGFfp,"%lf\t%lf\n",MSSpec[i].dMass,MSSpec[i].dIntensity);
		}
	}		
	fprintf(MGFfp,"END IONS\n");
	return true;
}

bool FTRawFile::XICCalibrateExt(double C[pdNum],FILE *MGFfp,long ScanNum)
{
	if(!GetMS2pMass(ScanNum)) return false;		
    //refine the parent mz value in the MS scan	
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return false;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;
	while(!PreFindMIntFtISO(it))//refine the prefine mass 
	{
		tempSNo--;
		if(tempSNo<SNFirst) return false;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);	
	}
	if(!GetFTingData(tempSNo,it)) return false;
	//it.pTIC=GetMSTIC();	
	//RTFromScanNum(tempSNo,&(it.RT));
	XICList_F.push_back(it);
	int Interupt_times=0;
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);			
		if(FindMIntFtISO(it)&&GetFTingData(tempSNo,it)) 
		{
			//it.pTIC=GetMSTIC();
			//RTFromScanNum(tempSNo,&(it.RT));
			XICList_F.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
	XICList.push_back(XICList_F[0]);
	Interupt_times=0;
	while(tempSNo<SNLast)
	{
		tempSNo++;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);		
		if(FindMIntFtISO(it)&&GetFTingData(tempSNo,it)) 
		{
			//it.pTIC=GetMSTIC();
			//RTFromScanNum(tempSNo,&(it.RT));
			XICList.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;	
		if(Interupt_times>=IT_LIM) break;
	}

	if(XICList.size()<=0) return false;

	vector<double> dMass;
	GetXICDataExt(&XICList,&dMass,C);
	double m,s;
	CalMST(&dMass,&m,&s);
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d_%lf.%d.%d\n",sFName,dMass.size(),s,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",m*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",m);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);
	
	if(ReadMS(ScanNum,MSSpec))
	{
		size=MSSpec.size();
		for(int i=0;i<size;i++)
		{
			if(MSSpec[i].dIntensity<1e-2) continue;
			fprintf(MGFfp,"%lf\t%lf\n",MSSpec[i].dMass,MSSpec[i].dIntensity);
		}
	}		
	fprintf(MGFfp,"END IONS\n");
	return true;
}

double FTRawFile::XICCalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI)
{		
    if(MSLevel(ScanNum)!=2) return 0;	
	//refine the parent mz value in the MS scan
	bMass=mzI;
	charge=ch;
	CalDefIsoDis(mzI*ch,IsoDisT);
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)==1) break;	
	}
	if(tempSNo<SNFirst) return ERROR_MASS;	
	if(!ReadMS(tempSNo,MSSpec)) return ERROR_MASS;	
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;	
	////////////////////////////
	if(!FindMIntFtISO(it)) return ERROR_MASS;
	if(!GetFTingData(tempSNo,it)) return ERROR_MASS;
	it.pTIC=GetMSTIC();	
	//it.ScanNum=1;//to record is MS/MS precusor
	it.ScanNum=tempSNo;
	it.charge=ch;

	RTFromScanNum(tempSNo,&(it.RT));
	XICList_F.push_back(it);
	//it.ScanNum=0;//to record is not MS/MS precusor
	int Interupt_times=0;
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);		
		if(FindMIntFtISO(it)&&GetFTingData(tempSNo,it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			it.ScanNum=tempSNo;
			if(it.goodness>0.9)	XICList_F.push_back(it);
			Interupt_times=0;			
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
    XICList.push_back(XICList_F[0]);
	Interupt_times=0;
	while(tempSNo<SNLast)
	{
		tempSNo++;
		if(MSLevel(tempSNo)!=1) continue;	
		ReadMS(tempSNo,MSSpec);				
		if(FindMIntFtISO(it)&&GetFTingData(tempSNo,it)) 
		{
			it.pTIC=GetMSTIC();
			RTFromScanNum(tempSNo,&(it.RT));
			it.ScanNum=tempSNo;
			if(it.goodness>0.9)	XICList.push_back(it);
			Interupt_times=0;	
		}
		else Interupt_times++;		
		if(Interupt_times>=IT_LIM) break;
	}

	vector<double> dMass;

	GetXICDataExt(&XICList,&dMass,C);
	double m,s;
	rCalMST(&dMass,&m,&s);	
	return m;
}


int FTRawFile::Calibrate(double C[5],char *MGFfile)
{
	FILE *MGFfp;
	MGFfp=fopen(MGFfile,"w");
	if(MGFfp==NULL) return 0;
	int CT=0;
	for(long scannum=SNFirst;scannum<=SNLast;scannum++)
	{
		if(!XICCalibrate(C,MGFfp,scannum))
		{
			CT+=ECalibrate(C,MGFfp,scannum);
		}
		else CT++;
	}
	fclose(MGFfp);
	return CT;
}

int FTRawFile::CalibrateExt(double C[pdNum],char *MGFfile)
{
	FILE *MGFfp;
	MGFfp=fopen(MGFfile,"w");
	if(MGFfp==NULL) return 0;
	int CT=0;
	for(long scannum=SNFirst;scannum<=SNLast;scannum++)
	{
		if(!XICCalibrateExt(C,MGFfp,scannum))
		{
			CT+=ECalibrateExt(C,MGFfp,scannum);
		}
		else CT++;
	}
	fclose(MGFfp);
	return CT;
}

double FTRawFile::CalibrateExt(double C[pdNum],long ScanNum,double mzinitial,int ch)
{
	double bReturn=XICCalibrateExt(C,ScanNum,ch,mzinitial);
	if(bReturn<0)
	{
		bReturn=ECalibrateExt(C,ScanNum,ch,mzinitial);
	}
	return bReturn;
}

int FTRawFile::SimpleCalibrate(double C[5],char *MGFfile)
{
	FILE *MGFfp;
	MGFfp=fopen(MGFfile,"w");
	if(MGFfp==NULL) return 0;
	int CT=0;
	for(long scannum=SNFirst;scannum<=SNLast;scannum++)
	{
		if(!SCalibrate(C,MGFfp,scannum))
		{
			CT+=ECalibrate(C,MGFfp,scannum);
		}
		else CT++;
	}
	fclose(MGFfp);
	return CT;
}

int FTRawFile::SimpleCalibrateExt(double C[pdNum],char *MGFfile)
{
	FILE *MGFfp;
	MGFfp=fopen(MGFfile,"w");
	if(MGFfp==NULL) return 0;
	int CT=0;
	for(long scannum=SNFirst;scannum<=SNLast;scannum++)
	{
		if(!SCalibrateExt(C,MGFfp,scannum))
		{
			CT+=ECalibrateExt(C,MGFfp,scannum);
		}
		else CT++;
	}
	fclose(MGFfp);
	return CT;
}

double FTRawFile::SimpleCalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI)
{
	double bR=SCalibrateExt(C,ScanNum,ch,mzI);
	if(bR<0)
	{
		bR=ECalibrateExt(C,ScanNum,ch,mzI);
	}
	return bR;
}

bool FTRawFile::GetTIC(long ScanNum,double *TIC, double *RT)
{
	*TIC=0;
	*RT=0;
	if(ScanNum>SNLast||ScanNum<SNFirst) return false;
	long tempSNo=ScanNum;	
	while(tempSNo>SNFirst)
	{
		tempSNo--;
		if(MSLevel(tempSNo)!=1) continue;	
		if(ReadMS(tempSNo,MSSpec))
		{
			size_t size=MSSpec.size();
			for(int i=0;i<size;i++)
			{
				*TIC+=MSSpec[i].dIntensity;
			}
			RTFromScanNum(tempSNo,RT);			
			break;
		}
	}
	return true;
}

void FTRawFile::CalDefIsoDis(double mass,double IsoAreaT[6])
{
	IsoAreaT[0]=0.9937*exp(-0.0006693 *mass);

	double item=0.0006956*mass-0.03438;
	IsoAreaT[1]=item*exp(-item);

	item=0.0006899 *mass-0.018;
	IsoAreaT[2]=item*item*exp(-item)/2;

	item=0.0006673*mass+0.03829;
	IsoAreaT[3]=item*item*item*exp(-item)/6;
	
	item=0.0006462*mass+0.1133;
	IsoAreaT[4]=item*item*item*item*exp(-item)/24;

	item=0.0006522*mass+0.1314;
	IsoAreaT[5]=item*item*item*item*item*exp(-item)/120;
}

void FTRawFile::EnableOutput(char *fname)
{
	fp_cal_data=fopen(fname,"w");
	if(fp_cal_data!=NULL) IsOutput=true;
	else IsOutput=false;
}

void FTRawFile::CloseOutput()
{
	if(IsOutput)
	{
		if(fp_cal_data!=NULL) fclose(fp_cal_data);
		IsOutput=false;
	}
}


double FTRawFile::PreFindMIntWt()
{	
	double dm=0.05;//the float point data reduction error range	
	size_t size=MSSpec.size();
	int middle=TFindL(MSSpec,0,size,bMass);
	if(middle>=size-1||middle<=0) return 0;

	int i=middle;
	//an increase adge, go up and down
	//                  /\
	//                 /  \
	//                /|   \
	//               / |    \
	//----------------------------------
	while(i<size-1)
	{
		//if(MSSpec[i+1].dIntensity<=1e-4) break;	
		if(MSSpec[i+1].dIntensity<=MSSpec[i].dIntensity)  break;
		if(MSSpec[i].dMass>dm+bMass) break;
		i++;
	}

	while(i<size-1)
	{
		//if(MSSpec[i+1].dIntensity<=1e-4) break;	
		if(MSSpec[i+1].dIntensity>MSSpec[i].dIntensity)  break;
		if(MSSpec[i].dMass>dm+bMass) break;
		i++;
	}
	//                  /\
	//                 /  \
	//                /   |\
	//               /    | \
	//----------------------------------
    //a decrease edge, go down directly
	int End=i;

	i=middle;
	while(i>0)
	{
		//if(MSSpec[i-1].dIntensity<=1e-4) break;	
		if(MSSpec[i].dIntensity>=MSSpec[i-1].dIntensity)  break;
		if(MSSpec[i].dMass+dm<bMass) break;
		i--;
	}

	while(i>0)
	{
		//if(MSSpec[i-1].dIntensity<=1e-4) break;	
		if(MSSpec[i].dIntensity<MSSpec[i-1].dIntensity)  break;
		if(MSSpec[i].dMass+dm<bMass) break;
		i--;
	}
	//int Begin=i;
	double massC=0;
	double sumI=0;
	for(/*i=Begin*/;i<=End;i++)
	{
		massC+=MSSpec[i].dIntensity*MSSpec[i].dMass;
		sumI+=MSSpec[i].dIntensity;
	}
	return massC/sumI;
}


double FTRawFile::PreFindMIntMP()
{	
	double dm=0.05;//the float point data reduction error range	
	size_t size=MSSpec.size();
	int middle=TFindL(MSSpec,0,size,bMass);
	if(middle>=size-1||middle<=0) return 0;

	int i=middle;
	//an increase adge, go up and down
	//                  /\
	//                 /  \
	//                /|   \
	//               / |    \
	//----------------------------------
	while(i<size-1)
	{
		//if(MSSpec[i+1].dIntensity<=1e-4) break;	
		if(MSSpec[i+1].dIntensity<=MSSpec[i].dIntensity)  break;
		if(MSSpec[i].dMass>dm+bMass) break;
		i++;
	}

	while(i<size-1)
	{
		//if(MSSpec[i+1].dIntensity<=1e-4) break;	
		if(MSSpec[i+1].dIntensity>MSSpec[i].dIntensity)  break;
		if(MSSpec[i].dMass>dm+bMass) break;
		i++;
	}
	//                  /\
	//                 /  \
	//                /   |\
	//               /    | \
	//----------------------------------
    //a decrease edge, go down directly
	int End=i;

	i=middle;
	while(i>0)
	{
		//if(MSSpec[i-1].dIntensity<=1e-4) break;	
		if(MSSpec[i].dIntensity>=MSSpec[i-1].dIntensity)  break;
		if(MSSpec[i].dMass+dm<bMass) break;
		i--;
	}

	while(i>0)
	{
		//if(MSSpec[i-1].dIntensity<=1e-4) break;	
		if(MSSpec[i].dIntensity<MSSpec[i-1].dIntensity)  break;
		if(MSSpec[i].dMass+dm<bMass) break;
		i--;
	}
	//int Begin=i;
	double massC=bMass;
	double maxI=0;
	for(/*i=Begin*/;i<=End;i++)
	{
		if(maxI<MSSpec[i].dIntensity)
		{
			massC=MSSpec[i].dMass;
			maxI=MSSpec[i].dIntensity;
		}
	}
	return massC;
}

double FTRawFile::GetSampleDm(int type,double mass)
{
	if(type==FT) return 6.06e-9*mass*mass-7.063e-11*mass+4.085e-8;//for FT
	else if(type==ORBITRAP) return 1.12e-009*mass*mass+2.204e-006*mass-0.0003366;//for Orbitrap
	else return 0;
}


#define ISO_W 2.0
//do not taken into account the pmass lost case
//that case was processed in the first finding as a special case
bool FTRawFile::PreFindMIntFtISO(IsoCluster &it)
{
	mPeak mt;
	if(!PreFindMIntFtNew()) return false;//determine bMass	
	if(!PreRefinePISO(it)) return false;//charge and bMass redetermined
	//ite is old program, not used here
	CalDefIsoDis(bMass*charge,IsoDisT);	//recalculation the iso dis for the final charge and bMass
	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		it.ISOP[i].dIntensity=0;
		it.ISOP[i].dMass=0;
		it.ISOP[i].rInt=0;
	}
	double backbmass=bMass;
	//back the bMass	;
	if(FindMIntFt(it.ISOP[0])==0)
	{
		bMass=backbmass;
		return false;
	}
	for(i=1;i<MAX_ISO;i++)
	{
		bMass=backbmass+ISO_DIFF[i]/charge;
		if(FindMIntFt(it.ISOP[i]) )break;
		if(it.ISOP[i].dIntensity<1e-6) break;
	}
	GetCOS(it,i);		
	bMass=backbmass;
	if(it.goodness<0.9) return false;
	return true;
}


//do not taken into account the pmass lost case
//that case was processed in the first finding as a special case
bool FTRawFile::FindMIntFtISO(IsoCluster &it)
{	
	if(charge<=0) return false;	
	if(IsoDisT[0]<=1e-4) return false;

	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		it.ISOP[i].dIntensity=0;
		it.ISOP[i].dMass=0;
		it.ISOP[i].rInt=0;
	}
	double backbmass=bMass;
	//back the bMass	;
	if(FindMIntFt(it.ISOP[0])==0)
	{
		bMass=backbmass;
		return false;
	}
	for(i=1;i<MAX_ISO;i++)
	{
		bMass=backbmass+ISO_DIFF[i]/charge;
		if(!FindMIntFt(it.ISOP[i])) break;
		if(it.ISOP[i].dIntensity<1e-6) break;
	}

	GetCOS(it,i);	
	bMass=backbmass;
	return true;
}

double FTRawFile::PreFindMIntFt()
{	
	//located to peak
	double dm=GetSampleDm(FT,bMass)*1.5;//the float point data reduction error range	
	int n=MSSpec.size();		
	int i,middle=LocatePmass(bMass-0.05,bMass+0.05,bMass);
	if(middle==-1)
	{	
		//externd the range
		middle=LocatePmassH(bMass-ISO_W,bMass+ISO_W,bMass);
		if(middle==-1) return bMass; //cna not find, use the original			
	}

	if(middle>=n) middle=n-1;
	mPeak ptold;
	ptold=MSSpec[middle];
	for(i=middle+1;i<n;i++)
	{			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass>ptold.dMass+dm)  break;	
		ptold=MSSpec[i];
	}	
	int End=i;//注意下面是小于，这里不需要减一

	i=middle;
	ptold=MSSpec[middle];
	while(i>0)
	{
		i--;			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass+dm<ptold.dMass)  break;	
		ptold=MSSpec[i];
	}
	i++;
	//int Begin=i;
	vector<double> pkmass;
	vector<double> pkintensity;
	vector<double> beta;
	//pkmass.push_back(LowR);
	//pkintensity.push_back(0);
	for(;i<End;i++)
	{		
		pkmass.push_back(MSSpec[i].dMass);
		pkintensity.push_back(MSSpec[i].dIntensity);
	}
	//pkmass.push_back(HighR);
	//pkintensity.push_back(0);
	double RSquare=Gauss_Fit(pkmass,pkintensity,beta);
	if(RSquare<0.9) return bMass;
	n=beta.size()/4;

	//find the most neibour peaks
	dm=0.1;	
	for(i=0;i<n;i++)
	{
		double tdm=abs(beta[4*i+1]-bMass);
		if(tdm<dm)
		{
			RSquare=beta[4*i+1];
			dm=tdm;
		}
	}
	return RSquare;	
}

bool FTRawFile::PreFindMIntFtNew()
{	
	//located to peak
	double dm=GetSampleDm(FT,bMass)*1.5;//the float point data reduction error range	
	int n=MSSpec.size();		
	int i,middle=LocatePmassH(bMass-0.05,bMass+0.05,bMass);
	if(middle==-1)
	{	
		//externd the range
		middle=LocatePmassH(bMass-ISO_W,bMass+ISO_W,bMass);
		if(middle==-1) return false; //cna not find, use the original			
	}

	if(middle>=n) middle=n-1;
	mPeak ptold;
	ptold=MSSpec[middle];
	for(i=middle+1;i<n;i++)
	{			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass>ptold.dMass+dm)  break;	
		ptold=MSSpec[i];
	}	
	int End=i;//注意下面是小于，这里不需要减一

	i=middle;
	ptold=MSSpec[middle];
	while(i>0)
	{
		i--;			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass+dm<ptold.dMass)  break;	
		ptold=MSSpec[i];
	}
	i++;
	vector<double> pkmass;
	vector<double> pkintensity;
	vector<double> beta;
	for(;i<End;i++)
	{		
		pkmass.push_back(MSSpec[i].dMass);
		pkintensity.push_back(MSSpec[i].dIntensity);
	}
	double RSquare=Gauss_Fit(pkmass,pkintensity,beta);
	if(RSquare<0.8) return false;
	n=beta.size()/4;

	//find the most neibour peaks
	double MaxInt=0;	
	for(i=0;i<n;i++)
	{
		if(MaxInt<beta[4*i])
		{	
			bMass=beta[4*i+1];
			MaxInt=beta[4*i];
		}	
	}
	return true;	
}


//bool FTRawFile::FindMass(double massE[5],long ScanNum,double massP,int ch)
//{
//	bMass=massP;
//	charge=ch;
//    //refine the parent mz value in the MS scan
//	long tempSNo=ScanNum;	
//	//seek to the last MS scan
//	while(tempSNo>SNFirst)
//	{
//		tempSNo--;
//		if(MSLevel(tempSNo)==1) break;	
//	}
//	if(tempSNo<SNFirst) return false;
//	if(!ReadMS(tempSNo,MSSpec)) return false;
//	//refine the pmz	
//	double bmass=PreFindMIntSim();
//	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d\n",sFName,ScanNum,charge);
//	//fprintf(MGFfp,"PMASS=%lf\n",cmass*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
//	fprintf(MGFfp,"PMASS=%lf\n",bmass);//maybe wrong
//	fprintf(MGFfp,"CHARGE=%d\n",charge);	
//	return true;
//}


int FTRawFile::LocatePmass(double begin,double end,double pmass)
{
	int i,middle=-1;
	double DM=0.02;
	int n=MSSpec.size();

	int iB=TFindL(0,n-1,begin);
	int iE=TFindH(0,n-1,end);

	for(i=iB;i<iE;i++)
	{
		if(MSSpec[i].dIntensity<=1e-6) continue;
		double Int1=0;
		double Int2=0;
		if(i>0) Int1=MSSpec[i-1].dIntensity;
		if(i<n-1) Int2=MSSpec[i+1].dIntensity;
		if(MSSpec[i].dIntensity>=Int1&&MSSpec[i].dIntensity>=Int2)
		{
			double dm=abs(MSSpec[i].dMass-pmass);
			if(dm<DM)
			{
				middle=i;
				DM=dm;
			}
		}	
	}
	return middle;
}

int FTRawFile::LocatePmass(double begin,double end,double pmass,double DM)
{
	int i,middle=-1;
	int n=MSSpec.size();

	int iB=TFindL(0,n-1,begin);
	int iE=TFindH(0,n-1,end);

	for(i=iB;i<iE;i++)
	{
		if(MSSpec[i].dIntensity<=1e-6) continue;
		double Int1=0;
		double Int2=0;
		if(i>0) Int1=MSSpec[i-1].dIntensity;
		if(i<n-1) Int2=MSSpec[i+1].dIntensity;
		if(MSSpec[i].dIntensity>=Int1&&MSSpec[i].dIntensity>=Int2)
		{
			double dm=abs(MSSpec[i].dMass-pmass);
			if(dm<DM)
			{
				middle=i;
				DM=dm;
			}
		}	
	}
	return middle;
}

int FTRawFile::LocatePmassH(double begin,double end,double pmass)
{
	int i,middle=-1;	
	int n=MSSpec.size();	
	double MaxInt=0;

	int iB=TFindL(0,n-1,begin);
	int iE=TFindH(0,n-1,end);

	for(i=iB;i<iE;i++)
	{		
		if(MSSpec[i].dIntensity<=1e-6) continue;		
		if(MaxInt<MSSpec[i].dIntensity)
		{
			//pmass=MSSpec[i].dMass;
			MaxInt=MSSpec[i].dIntensity;
			middle=i;
		}
	}
	return middle;
}

//the user must make sure the index no t exceed the band
int FTRawFile::TFindL(int b,int e,double fmz)
{
	if(e-b<=1) return b;
	if(MSSpec[b].dMass<=fmz&&MSSpec[e].dMass>=fmz) return b;
	int idx=(b+e)/2;	
	if(MSSpec[idx].dMass>fmz)
	{
		return TFindL(b,idx,fmz);
	}
	else return TFindL(idx,e,fmz);
}

int FTRawFile::TFindH(int b,int e,double fmz)
{
	if(e-b<=1) return e;
	if(MSSpec[b].dMass<=fmz&&MSSpec[e].dMass>=fmz) return e;
	int idx=(b+e)/2;
	if(MSSpec[idx].dMass>fmz)
	{
		return TFindH(b,idx,fmz);
	}
	else return TFindH(idx,e,fmz);
}


int FTRawFile::PreFindMIntFt(mPeak &mt)
{	
	//located to peak	
	double dm=GetSampleDm(FT,bMass);//the float point data reduction error range	
	double Pre_WD=dm*10;
	// for the Rounding error
	if(Pre_WD<0.05) Pre_WD=0.05;
	dm*=1.5;
	int i,middle=LocatePmass(bMass-Pre_WD,bMass+Pre_WD,bMass,Pre_WD);
	int bReturn=0;	
	if(middle==-1)
	{	
		//externd the range
		middle=LocatePmassH(bMass-ISO_W,bMass+ISO_W,bMass);
		if(middle==-1) return 0; //can not find, use the original	
		bReturn|=0x08;
	}
	bReturn|=0x01;

	size_t n=MSSpec.size();
	if(middle>=n) middle=n-1;
	int idx_old=middle;
	for(i=middle+1;i<n;i++)
	{			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass>MSSpec[idx_old].dMass+dm)  break;	
		idx_old=i;
	}	
	int End=i;//注意下面是小于，这里不需要减一

	i=middle;
	idx_old=middle;
	while(i>0)
	{
		i--;			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass+dm<MSSpec[idx_old].dMass)  break;	
		idx_old=i;
	}
	i++;
	//int Begin=i;
	vector<double> pkmass;
	vector<double> pkintensity;
	vector<double> beta;
	//pkmass.push_back(LowR);
	//pkintensity.push_back(0);
	for(;i<End;i++)
	{		
		pkmass.push_back(MSSpec[i].dMass);
		pkintensity.push_back(MSSpec[i].dIntensity);
	}
	//pkmass.push_back(HighR);
	//pkintensity.push_back(0);
	mt.RSquare=Gauss_Fit(pkmass,pkintensity,beta);
	if(mt.RSquare<0.9) bReturn|=0x02;
	n=beta.size()/4;

	//find the most neibour peaks
	dm=0.1;	
	double MaxInt=0;
	int idx=-1;
	for(i=0;i<n;i++)
	{
		double tdm=abs(beta[4*i+1]-bMass);
		if(MaxInt<beta[4*i]) MaxInt=beta[4*i];
		if(tdm<dm)
		{
			idx=i;		
			dm=tdm;
		}
	}
	if(idx!=-1)
	{
		mt.AssignFit(beta,4*idx);	
		//added on 2010.11.10
		double bt[4];
		for(i=0;i<4;i++) bt[i]=beta[idx*4+i];
		CalSKPar(pkmass,pkintensity,bt,mt.skew);
		//
	}
	n=MSSpec.size();
	MaxInt=0;
	for(i=0;i<n;i++)
	{
		if(MaxInt<MSSpec[i].dIntensity) MaxInt=MSSpec[i].dIntensity;
	}
	mt.rInt=mt.dIntensity/MaxInt;
	if(mt.dIntensity<MaxInt) bReturn|=0x04;	
	return bReturn;	
}

//only work for those with one peaks, beta size is 4
////four parameters sum of +/- before and after the peak
	//                  /|\
	//             +/- / | \+/-
	//                /  |  \
	//               /   |   \
	//----------------------------------
//x.size must be equal to y.size, this routine doese not check
void FTRawFile::CalSKPar(vector<double> x,vector<double>y,double beta[4],double ReturnPar[4])
{
	size_t n=x.size();
    size_t i;
	ReturnPar[0]=0;
	ReturnPar[1]=0;
	ReturnPar[2]=0;
	ReturnPar[3]=0;
	for(i=0;i<n;i++)
	{
		double ytemp=Gauss_F(x[i],beta);
		//double sumy=(ytemp+y[i])/2;
		double sumy=y[i];//modified on 2010.12.16
		if(sumy<=1e-4) continue;
		ytemp=y[i]-ytemp;
		ytemp/=sumy;
		if(ytemp<0)
		{
			if(x[i]>beta[1]) ReturnPar[0]+=ytemp;
			else ReturnPar[2]+=ytemp;
		}
		else 
		{
			if(x[i]>beta[1]) ReturnPar[1]+=ytemp;
			else ReturnPar[3]+=ytemp;
		}
	}
}

void FTRawFile::CalSKPar(vector<double> x,vector<double>y,vector<double> beta,double ReturnPar[4])
{
	size_t n=x.size();
    size_t i;
	ReturnPar[0]=0;
	ReturnPar[1]=0;
	ReturnPar[2]=0;
	ReturnPar[3]=0;
	for(i=0;i<n;i++)
	{
		double ytemp=Gauss_F(x[i],beta);
		//double sumy=(ytemp+y[i])/2;
		double sumy=y[i];//modified on 2010.12.16
		if(sumy<=1e-4) continue;
		ytemp-=y[i];
		ytemp/=sumy;
		if(ytemp<0)
		{
			if(x[i]>beta[1]) ReturnPar[0]+=ytemp;
			else ReturnPar[2]+=ytemp;
		}
		else 
		{
			if(x[i]>beta[1]) ReturnPar[1]+=ytemp;
			else ReturnPar[3]+=ytemp;
		}
	}
}

int FTRawFile::FindMIntFt(mPeak &mt)
{	
	//located to peak	
	double dm=GetSampleDm(FT,bMass);//the float point data reduction error range	
	double Pre_WD=dm*10;
	dm*=1.5;
	int i,middle=LocatePmass(bMass-Pre_WD,bMass+Pre_WD,bMass,Pre_WD);
	int bReturn=0;	
	if(middle==-1) return 0; //can not find, use the original	
	bReturn|=0x01;
	size_t n=MSSpec.size();
	if(middle>=n) middle=n-1;
	int idx_old=middle;
	for(i=middle+1;i<n;i++)
	{			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass>MSSpec[idx_old].dMass+dm)  break;	
		idx_old=i;
	}	
	int End=i;//注意下面是小于，这里不需要减一

	i=middle;
	idx_old=middle;
	while(i>0)
	{
		i--;			
		if(MSSpec[i].dIntensity<=1e-4) break;	
		if(MSSpec[i].dMass+dm<MSSpec[idx_old].dMass)  break;	
		idx_old=i;
	}
	i++;
	//int Begin=i;
	vector<double> pkmass;
	vector<double> pkintensity;
	vector<double> beta;
	//pkmass.push_back(LowR);
	//pkintensity.push_back(0);
	for(;i<End;i++)
	{		
		pkmass.push_back(MSSpec[i].dMass);
		pkintensity.push_back(MSSpec[i].dIntensity);
	}
	//pkmass.push_back(HighR);
	//pkintensity.push_back(0);
	mt.RSquare=Gauss_Fit(pkmass,pkintensity,beta);
	if(mt.RSquare<0.9) bReturn|=0x02;
	n=beta.size()/4;

	//find the most neibour peaks
	dm=0.1;	
	double MaxInt=0;
	int idx=-1;
	for(i=0;i<n;i++)
	{
		double tdm=abs(beta[4*i+1]-bMass);
		if(MaxInt<beta[4*i]) MaxInt=beta[4*i];
		if(tdm<dm)
		{
			idx=i;
			//mt.AssignFit(beta,4*i);			
			dm=tdm;
		}
	}
	if(idx!=-1)
	{
		mt.AssignFit(beta,4*idx);	
		//added on 2010.11.10
		double bt[4];
		for(i=0;i<4;i++) bt[i]=beta[idx*4+i];
		CalSKPar(pkmass,pkintensity,bt,mt.skew);
		//
	}
	if(mt.dIntensity<MaxInt) bReturn|=0x04;	
	n=MSSpec.size();
	MaxInt=0;
	for(i=0;i<n;i++)
	{
		if(MaxInt<MSSpec[i].dIntensity) MaxInt=MSSpec[i].dIntensity;
	}
	mt.rInt=mt.dIntensity/MaxInt;

	return bReturn;	
}


void FTRawFile::GetCOS(IsoCluster &it,int num)
{	
	int i,numMatch=0;	
	for(i=0;i<num;i++)
	{
		if(it.ISOP[i].dIntensity>1e-6)numMatch++;
	}	
	double Intsum=0;
	double Intsum1=0;
	double Intsum2=0;
	for(i=0;i<num;i++)
	{
		Intsum+=it.ISOP[i].dIntensity*it.ISOP[i].dIntensity;
		Intsum1+=IsoDisT[i]*IsoDisT[i];
		Intsum2+=it.ISOP[i].dIntensity*IsoDisT[i];
	}
	Intsum=sqrt(Intsum1*Intsum);
	if(Intsum<=1e-6) it.goodness=0;	
	else it.goodness=Intsum2/Intsum;
}

double FTRawFile::BasePeakInt()
{
	int s=MSSpec.size();
	double bs=0;
	for(int i=0;i<s;i++)
	{
		if(bs<MSSpec[i].dIntensity) bs=MSSpec[i].dIntensity;
	}
	return bs;
}
//added on 2010.10.13
bool FTRawFile::GetFTingData(double pMass,long ScanNum,CalData &dt,int ch)
{
	bMass=pMass;
	charge=ch;
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		if(MSLevel(tempSNo)==1) break;
		tempSNo--;	
	}
	if(tempSNo<SNFirst) return 0;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	dt.tic=GetMSTIC();	
	RTFromScanNum(tempSNo,&(dt.RT));	
	IsoCluster it;	
	////////////////////////////
	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		IsoDisT[i]=dt.IsotopicT[i];
	}
	if(!FindMIntFtISO(it)) return false;
	dt.mzeF=dt.mze;//removealbe fro formal use, added by zhangjy on 2010.11.2
	dt.mze=it.ISOP[0].dMass;
	dt.mzeXC=GetPMassL(tempSNo,dt.mzeF);//removealbe fro formal use, added by zhangjy on 2010.11.2
	//added on 2010.12.16
	for(i=0;i<4;i++) dt.PShape[i]=it.ISOP[0].skew[i];
	//////////////////////////////////compare infor
	for(i=0;i<MAX_ISO;i++)       
	{
		dt.IsotopicE[i]=it.ISOP[i].dIntensity;
	}
	dt.aInt=it.ISOP[0].dIntensity;
	dt.rInt=it.ISOP[0].rInt;	
	//added on 2010.10.11
	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusData(tempSNo,dt);
}

bool FTRawFile::GetFTingDataAll(double pMass,long ScanNum,CalData &dt,int ch)
{
	bMass=pMass;
	charge=ch;
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		if(MSLevel(tempSNo)==1) break;
		tempSNo--;	
	}
	if(tempSNo<SNFirst) return 0;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	dt.tic=GetMSTIC();	
	RTFromScanNum(tempSNo,&(dt.RT));	
	IsoCluster it;	
	////////////////////////////
	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		IsoDisT[i]=dt.IsotopicT[i];
	}
	if(!FindMIntFtISO(it)) return false;
	dt.mzeF=dt.mze;//removealbe fro formal use, added by zhangjy on 2010.11.2
	dt.mze=it.ISOP[0].dMass;
	dt.mzeXC=GetPMassL(tempSNo,dt.mzeF);//removealbe fro formal use, added by zhangjy on 2010.11.2
	//added on 2010.12.16
	for(i=0;i<4;i++) dt.PShape[i]=it.ISOP[0].skew[i];
	//////////////////////////////////compare infor
	for(i=0;i<MAX_ISO;i++)       
	{
		dt.IsotopicE[i]=it.ISOP[i].dIntensity;
	}
	dt.aInt=it.ISOP[0].dIntensity;
	dt.rInt=it.ISOP[0].rInt;	
	//added on 2010.10.11
	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusDataAll(tempSNo,dt);
}

bool FTRawFile::GetFTingData(double pMass,long ScanNum,CalData &dt)
{
	bMass=pMass;
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{		
		if(MSLevel(tempSNo)==1) break;	
		tempSNo--;
	}
	if(tempSNo<SNFirst) return 0;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	dt.tic=GetMSTIC();	
	RTFromScanNum(tempSNo,&(dt.RT));
	mPeak it;
	if(!FindMIntFt(it)) return false;
	dt.mze=it.dMass;
	//added on 2010.12.16
	for(int i=0;i<4;i++) dt.PShape[i]=it.skew[i];
	//added on 2010.10.11
	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusData(tempSNo,dt);
	//return true;
}

bool FTRawFile::GetFTingDataAll(double pMass,long ScanNum,CalData &dt)
{
	bMass=pMass;
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{		
		if(MSLevel(tempSNo)==1) break;	
		tempSNo--;
	}
	if(tempSNo<SNFirst) return 0;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	dt.tic=GetMSTIC();	
	RTFromScanNum(tempSNo,&(dt.RT));
	mPeak it;
	if(!FindMIntFt(it)) return false;
	dt.mze=it.dMass;
	//added on 2010.12.16
	for(int i=0;i<4;i++) dt.PShape[i]=it.skew[i];
	//added on 2010.10.11
	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusDataAll(tempSNo,dt);
	//return true;
}

bool FTRawFile::GetFTingData(long ScanNum,CalData &dt)
{	
    //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		if(MSLevel(tempSNo)==1) break;	
		tempSNo--;
	}
	if(tempSNo<SNFirst) return 0;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	dt.tic=GetMSTIC();	
	RTFromScanNum(tempSNo,&(dt.RT));
	mPeak it;
	if(!FindMIntFt(it)) return false;
	dt.mze=it.dMass;	
	//added on 2010.12.16
	for(int i=0;i<4;i++) dt.PShape[i]=it.skew[i];
	//added on 2010.10.11
	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusData(tempSNo,dt);
	//return true;
}

bool FTRawFile::GetFTingData(double pMass,long ScanNum,IsoCluster &dt)
{
	bMass=pMass;
	 //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{		
		if(MSLevel(tempSNo)==1) break;	
		tempSNo--;
	}
	if(tempSNo<SNFirst) return 0;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	dt.pTIC=GetMSTIC();	
	RTFromScanNum(tempSNo,&(dt.RT));	
	//added on 2010.10.11
	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusData(tempSNo,dt);
	//return true;
}

bool FTRawFile::GetFTingData(long ScanNum,IsoCluster &dt)
{
	 //refine the parent mz value in the MS scan
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{		
		if(MSLevel(tempSNo)==1) break;	
		tempSNo--;
	}
	if(tempSNo<SNFirst) return 0;
	if(!ReadMS(tempSNo,MSSpec)) return false;	
	//refine the pmz
	dt.pTIC=GetMSTIC();	
	RTFromScanNum(tempSNo,&(dt.RT));	
	//added on 2010.10.11
	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusData(tempSNo,dt);
	//return true;
}

int FTRawFile::GetParIdx(CString &sLabel)
{
	for(int i=0;i<MAX_PAR_NUM;i++)
	{
		if(sLabel.Find(Parname[i].c_str())==0) return i;
	}
	return -1;
}


int FTRawFile::GetParIdxAll(CString &sLabel)
{
	for(int i=0;i<MAX_PAR_NUM_ALL;i++)
	{
		if(sLabel.Find(ParnameALL[i].c_str())==0) return i;
	}
	return -1;
}

bool FTRawFile::GetFTStatusData(long ScanNum,CalData &dt)
{
	double dStatusLogRT = 0.0;
	VARIANT varLabels;
	VariantInit(&varLabels);
	VARIANT varValues;
	VariantInit(&varValues);
	long nArraySize = 0;

	GetStatusLogForScanNum(ScanNum, 
	&dStatusLogRT, 
	&varLabels, 
	&varValues, 
	&nArraySize);
	if(nArraySize<=0) return false;

	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psaLabels = varLabels.parray;
	varLabels.parray = NULL;

	SAFEARRAY FAR* psaValues = varValues.parray;
	varValues.parray = NULL;

	BSTR* pbstrLabels = NULL;
	BSTR* pbstrValues = NULL;

	if(FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels))))
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		//printf("Failed to access labels array");
		return false;
	}

	if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		SafeArrayUnaccessData(psaValues);
		SafeArrayDestroy(psaValues);
		//printf("Failed to access values array");
		return false;
	}

	CString sLabel;
	CString sData;
	long i;
	for(i=0;i<MAX_PAR_NUM;i++) dt.FTStatus[i]=-1000;

	bool isR=false;	
	for( i=0;i<nArraySize;i++ )
	{
		sLabel = pbstrLabels[i];
		sData = pbstrValues[i];
		int idx=GetParIdx(sLabel);
		if(idx!=-1) 
		{
			if(idx==Repeat_IDX)
			{
				if(isR) idx++;
				else isR=true;
			}
			dt.FTStatus[idx]=atof((LPCTSTR)sData);	
		}
		//printf("Name=%s\n",(LPCTSTR)sLabel);
		//printf("Value=%s\n",(LPCTSTR)sData);
	}
	// Delete the SafeArray
	SafeArrayUnaccessData(psaLabels);
	SafeArrayDestroy(psaLabels);
	SafeArrayUnaccessData(psaValues);
	SafeArrayDestroy(psaValues);	
	return true;
}

bool FTRawFile::GetFTStatusDataAll(long ScanNum,CalData &dt)
{
	double dStatusLogRT = 0.0;
	VARIANT varLabels;
	VariantInit(&varLabels);
	VARIANT varValues;
	VariantInit(&varValues);
	long nArraySize = 0;

	GetStatusLogForScanNum(ScanNum, 
	&dStatusLogRT, 
	&varLabels, 
	&varValues, 
	&nArraySize);
	if(nArraySize<=0) return false;

	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psaLabels = varLabels.parray;
	varLabels.parray = NULL;

	SAFEARRAY FAR* psaValues = varValues.parray;
	varValues.parray = NULL;

	BSTR* pbstrLabels = NULL;
	BSTR* pbstrValues = NULL;

	if(FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels))))
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		//printf("Failed to access labels array");
		return false;
	}

	if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		SafeArrayUnaccessData(psaValues);
		SafeArrayDestroy(psaValues);
		//printf("Failed to access values array");
		return false;
	}

	CString sLabel;
	CString sData;
	long i;
	for(i=0;i<MAX_PAR_NUM_ALL;i++) dt.FTStatus[i]=-1000;

	bool isR=false;	
	for( i=0;i<nArraySize;i++ )
	{
		sLabel = pbstrLabels[i];
		sData = pbstrValues[i];
		int idx=GetParIdxAll(sLabel);
		if(idx!=-1) 
		{
			if(idx==Repeat_IDX_GLOBAL)
			{
				if(isR) idx++;
				else isR=true;
			}
			dt.FTStatus[idx]=atof((LPCTSTR)sData);	
		}
		//printf("Name=%s\n",(LPCTSTR)sLabel);
		//printf("Value=%s\n",(LPCTSTR)sData);
	}
	// Delete the SafeArray
	SafeArrayUnaccessData(psaLabels);
	SafeArrayDestroy(psaLabels);
	SafeArrayUnaccessData(psaValues);
	SafeArrayDestroy(psaValues);	
	return true;
}

bool FTRawFile::GetFTEJTData(long ScanNum,CalData &dt)
{
	long nArraySize = 0;
	VARIANT varLabels;
	VariantInit(&varLabels);
	VARIANT varValues;
	VariantInit(&varValues);

	GetTrailerExtraForScanNum(ScanNum, 
	&varLabels, 
	&varValues, 
	&nArraySize);
	if(nArraySize<=0) return false;
	
	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psaLabels = varLabels.parray;
	varLabels.parray = NULL;

	SAFEARRAY FAR* psaValues = varValues.parray;
	varValues.parray = NULL;

	BSTR* pbstrLabels = NULL;
	BSTR* pbstrValues = NULL;

	if( FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels) ) ) )
	{
		SafeArrayUnaccessData( psaLabels );
		SafeArrayDestroy( psaLabels );
		//printf("Failed to access labels array");
		return false;
	}

	if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
	{
		SafeArrayUnaccessData( psaLabels );
		SafeArrayDestroy( psaLabels );
		SafeArrayUnaccessData( psaValues );
		SafeArrayDestroy( psaValues );
		//printf("Failed to access values array");
		return false;
	}

	CString sLabel;
	CString sData;
	for( long i=0; i<nArraySize; i++ )
	{
		sLabel = pbstrLabels[i];
		sData = pbstrValues[i];

		if(sLabel.Find("Ion Injection Time")!=-1)
		{
			dt.IIt=atof((LPCTSTR)sData);
		}
		else if(sLabel.Find("Elapsed Scan Time")!=-1)
		{
			dt.ESt=atof((LPCTSTR)sData);
		}		

		//printf("Name=%s\n",(LPCTSTR)sLabel);
		//printf("Value=%s\n",(LPCTSTR)sData);
	}
	//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",RT,IIT,EST,dBasePeakMass,dBasePeakIntensity,dTIC);

	// Delete the SafeArray
	SafeArrayUnaccessData(psaLabels);
	SafeArrayDestroy(psaLabels);
	SafeArrayUnaccessData(psaValues);
	SafeArrayDestroy(psaValues);	
	return true;
}

//ScanNum may be a MS2 index, will be convert the cossponding MS index automatically
bool FTRawFile::GetFTExtraData(long ScanNum,CalData &dt)
{
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{		
		if(MSLevel(tempSNo)==1) break;	
		tempSNo--;
	}
	if(tempSNo<SNFirst) return false;

	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusData(tempSNo,dt);
}

bool FTRawFile::GetFTEJTData(long ScanNum,IsoCluster &dt)
{
	long nArraySize = 0;
	VARIANT varLabels;
	VariantInit(&varLabels);
	VARIANT varValues;
	VariantInit(&varValues);

	GetTrailerExtraForScanNum(ScanNum, 
	&varLabels, 
	&varValues, 
	&nArraySize);
	if(nArraySize<=0) return false;
	
	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psaLabels = varLabels.parray;
	varLabels.parray = NULL;

	SAFEARRAY FAR* psaValues = varValues.parray;
	varValues.parray = NULL;

	BSTR* pbstrLabels = NULL;
	BSTR* pbstrValues = NULL;

	if( FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels) ) ) )
	{
		SafeArrayUnaccessData( psaLabels );
		SafeArrayDestroy( psaLabels );
		//printf("Failed to access labels array");
		return false;
	}

	if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
	{
		SafeArrayUnaccessData( psaLabels );
		SafeArrayDestroy( psaLabels );
		SafeArrayUnaccessData( psaValues );
		SafeArrayDestroy( psaValues );
		//printf("Failed to access values array");
		return false;
	}

	CString sLabel;
	CString sData;
	for( long i=0; i<nArraySize; i++ )
	{
		sLabel = pbstrLabels[i];
		sData = pbstrValues[i];

		if(sLabel.Find("Ion Injection Time")!=-1)
		{
			dt.IIt=atof((LPCTSTR)sData);
		}
		else if(sLabel.Find("Elapsed Scan Time")!=-1)
		{
			dt.ESt=atof((LPCTSTR)sData);
		}		

		//printf("Name=%s\n",(LPCTSTR)sLabel);
		//printf("Value=%s\n",(LPCTSTR)sData);
	}
	//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",RT,IIT,EST,dBasePeakMass,dBasePeakIntensity,dTIC);

	// Delete the SafeArray
	SafeArrayUnaccessData(psaLabels);
	SafeArrayDestroy(psaLabels);
	SafeArrayUnaccessData(psaValues);
	SafeArrayDestroy(psaValues);	
	return true;
}

//ScanNum may be a MS2 index, will be convert the cossponding MS index automatically
bool FTRawFile::GetFTExtraData(long ScanNum,IsoCluster &dt)
{
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>SNFirst)
	{
		if(MSLevel(tempSNo)==1) break;
		tempSNo--;	
	}
	if(tempSNo<SNFirst) return false;

	if(!GetFTEJTData(tempSNo,dt)) return false;
	return GetFTStatusData(tempSNo,dt);
}

bool FTRawFile::GetFTStatusData(long ScanNum,IsoCluster &dt)
{
	double dStatusLogRT = 0.0;
	VARIANT varLabels;
	VariantInit(&varLabels);
	VARIANT varValues;
	VariantInit(&varValues);
	long nArraySize = 0;

	GetStatusLogForScanNum(ScanNum, 
	&dStatusLogRT, 
	&varLabels, 
	&varValues, 
	&nArraySize);
	if(nArraySize<=0) return false;

	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psaLabels = varLabels.parray;
	varLabels.parray = NULL;

	SAFEARRAY FAR* psaValues = varValues.parray;
	varValues.parray = NULL;

	BSTR* pbstrLabels = NULL;
	BSTR* pbstrValues = NULL;

	if(FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels))))
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		//printf("Failed to access labels array");
		return false;
	}

	if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		SafeArrayUnaccessData(psaValues);
		SafeArrayDestroy(psaValues);
		//printf("Failed to access values array");
		return false;
	}

	CString sLabel;
	CString sData;
	long i;
	for(i=0;i<MAX_PAR_NUM;i++) dt.FTStatus[i]=-1000;

	bool isR=false;	
	for( i=0;i<nArraySize;i++ )
	{
		sLabel = pbstrLabels[i];
		sData = pbstrValues[i];
		int idx=GetParIdx(sLabel);
		if(idx!=-1) 
		{
			if(idx==Repeat_IDX)
			{
				if(isR) idx++;
				else isR=true;
			}
			dt.FTStatus[idx]=atof((LPCTSTR)sData);	
		}
		//printf("Name=%s\n",(LPCTSTR)sLabel);
		//printf("Value=%s\n",(LPCTSTR)sData);
	}
	// Delete the SafeArray
	SafeArrayUnaccessData(psaLabels);
	SafeArrayDestroy(psaLabels);
	SafeArrayUnaccessData(psaValues);
	SafeArrayDestroy(psaValues);	
	return true;
}
