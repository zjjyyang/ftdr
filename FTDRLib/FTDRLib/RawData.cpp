#include "stdafx.h"
#include "RawData.h"
#include "OutlierRem.h"
#include "OTrace.h"
#include "IntRModel.h"
#include "GaussFit.h"
#include "mzXMLOut.h"
#include "MassMatrixR.h"

MS2Header::MS2Header()
{
	pmz=0;
	scan=0;
	IW=0;
	MS1Scan=-1;

}

MS2Header::~MS2Header()
{

}

MS2Header::MS2Header(const MS2Header &cp)
{
	pmz=cp.pmz;
	scan=cp.scan;
	IW=cp.IW;
	ch=cp.ch;
	MS1Scan=cp.MS1Scan;
}

MS2Header &MS2Header::operator =(const MS2Header &cp)
{
	pmz=cp.pmz;
	scan=cp.scan;
	IW=cp.IW;
	ch=cp.ch;
	MS1Scan=cp.MS1Scan;
	return *this;
}

void MS2Header::Write2File(FILE *fp)
{
	fwrite(&pmz,sizeof(double),1,fp);
	fwrite(&scan,sizeof(long),1,fp);
	fwrite(&IW,sizeof(double),1,fp);
	fwrite(&ch,sizeof(int),1,fp);
	fwrite(&MS1Scan,sizeof(long),1,fp);
}

void MS2Header::ReadFromFile(FILE *fp)
{
	fread(&pmz,sizeof(double),1,fp);
	fread(&scan,sizeof(long),1,fp);
	fread(&IW,sizeof(double),1,fp);
	fread(&ch,sizeof(int),1,fp);
	fread(&MS1Scan,sizeof(long),1,fp);
}

CalReturn::CalReturn()
{
	cal_mean=0;
	cal_std=0;
	charge=0;
	rel_int=0;
	iso_goodness=0;
	iso_num=0;
	abs_int=0;
}

CalReturn::~CalReturn()
{

}

CalReturn::CalReturn(const CalReturn &cp)
{
	cal_mean=cp.cal_mean;
	cal_std=cp.cal_std;
	charge=cp.charge;
	rel_int=cp.rel_int;
	iso_goodness=cp.iso_goodness;
	iso_num=cp.iso_num;
	spectrum_num=cp.spectrum_num;
	abs_int=cp.abs_int;
}

CalReturn &CalReturn::operator=(const CalReturn &cp)
{
	cal_mean=cp.cal_mean;
	cal_std=cp.cal_std;
	charge=cp.charge;
	rel_int=cp.rel_int;
	iso_goodness=cp.iso_goodness;
	iso_num=cp.iso_num;
	spectrum_num=cp.spectrum_num;
	abs_int=cp.abs_int;
	return *this;
}

CalReturn &CalReturn::operator=(const CalData &cp)
{
	charge=cp.ch;
	iso_goodness=cp.goodness;
	iso_num=cp.IsoNum;
	if(cp.BaseInt>1e-6)	rel_int=cp.IsotopicE[0]/cp.BaseInt;
	else rel_int=0;
	abs_int=cp.GetmaxInt();
	return *this;
}

RawData::RawData(void)
{
/*	current_mz=0;
	current_ch=0;*/	
	raw_name[0]='\0';
	for(int i=0;i<MAX_ISO;i++) IsoDisT[i]=0;
	MET=15;
	CalModel=NULL;
	fp_cal_data=NULL;
	instrument_T=INSTRUMENTMODEL_UNDEF;
	TIF=NULL;
	Error_code=0;
	CalTotal=0;
	modeTotal=0;
	ppmLV=-10;
	ppmHV=10;

	METModel[0]=0;
	METModel[1]=0;
	METModel[2]=0;
	METModel[3]=0;
	IsMETModel=false;

	IntRModelData[0]='\0';
	OIntR_Data=false;
	IsOutput=false;

	ppbLevel=0;

	//	
	MzIntModel[0]=1.12e-009;
	MzIntModel[1]=2.204e-006;
	MzIntModel[2]=-0.0003366;
	MzIntModel[3]=1.0;

	IsUseExtendPar=true;

#ifdef CASE_TEST
	ParentIons=NULL;
#endif

	ModelType=MODEL_SVM;

	TailerIdx[0]=0;
	TailerIdx[1]=0;

	ms2format=MS2MGF;
}

RawData::~RawData(void)
{	
	if(IsOutput)
	{
		if(fp_cal_data!=NULL) fclose(fp_cal_data);		
	}
	if(TIF!=NULL) delete TIF;
}

long RawData::seekScan(long scannum,size_t a,size_t b)
{
	if(b<a) return -1;
	if(b-a<=1) return a;	
	int idx=(b+a)/2;
	if(ScanMap.at(idx)>=scannum)
	{
		return seekScan(scannum,a,idx);
	}
	else return seekScan(scannum,idx,b);
}


//adding version,mzt is assigned inside
bool RawData::PreCalData(CalData &ct,double mze,int ch,int MS2Scan,double mzt)
{	
	size_t ScanTotal=ScanMap.size();
	int MS1Scan=seekScan(MS2Scan,0,ScanTotal);
	int tempIDX=MS1Scan;	
	double RTOld=MSList[MS1Scan].RT;
	ct.mzt=mzt;
	if(tempIDX>=0)
	{
		if(MSList[tempIDX].FindISO(ct,ch,mze))	return true;		
	}
	return false;	
}

//adding version,mzt is assigned inside
void RawData::GetXICData(vector<CalData> &XICList,double mze,int ch,int MS2Scan,double mzt)
{	
	size_t ScanTotal=ScanMap.size();
	int MS1Scan=seekScan(MS2Scan,0,ScanTotal);
	int tempIDX=MS1Scan;
	CalData ct;
	double RTOld=MSList[MS1Scan].RT;
	ct.mzt=mzt;
	int interput=0;

	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO(ct,ch,mzt);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;		
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{
			XICList.push_back(ct);
			interput=0;
			RTOld=MSList[tempIDX].RT;			
		}
		
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;

	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO(ct,ch,mzt);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;	
			if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			XICList.push_back(ct);
			interput=0;
			RTOld=MSList[tempIDX].RT;		
		}		
		tempIDX++;		
	}
}

//adding version,mzt is assigned inside
void RawData::GetXICData_dMET(vector<CalData> &XICList,double mze,int ch,int MS2Scan,double mzt)
{	
	size_t ScanTotal=ScanMap.size();
	int MS1Scan=seekScan(MS2Scan,0,ScanTotal);
	int tempIDX=MS1Scan;
	CalData ct;
	double RTOld=MSList[MS1Scan].RT;
	ct.mzt=mzt;
	int interput=0;
	vector<CalData> XICListL;
	vector<CalData> XICListR;
	size_t i,TryT=RLEXTEND+LLEXTEND;
	i=0;

	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO_dMET(ct,ch,mzt,ppmLV,ppmHV);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{			
			XICListL.push_back(ct);		
			interput=0;
			RTOld=MSList[tempIDX].RT;
			i++;
			if(TryT>0&&i>=TryT) break;
		}		
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;
	interput=0;	
	RTOld=MSList[MS1Scan].RT;
	i=0;
	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO_dMET(ct,ch,mzt,ppmLV,ppmHV);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			XICListR.push_back(ct);
			interput=0;
			RTOld=MSList[tempIDX].RT;	
			i++;
			if(TryT>0&&i>=TryT) break;
		}
		
		tempIDX++;		
	}
	size_t countL=XICListL.size();
	size_t countR=XICListR.size();
	if(RLEXTEND>0&&LLEXTEND>0)
	{
		if(countR<RLEXTEND)  
		{
			i=TryT-countR;
			if(countL>i) countL=i;
		}
		else if(countL<LLEXTEND) 
		{
			i=TryT-countL;
			if(countR>i) countR=i;		
		}
		else
		{
			countR=RLEXTEND;
			countL=RLEXTEND;
		}
	}
	
	for(i=0;i<countL;i++)	XICList.push_back(XICListL[i]);	
	for(i=0;i<countR;i++)	XICList.push_back(XICListR[i]);	
	//printf("Left add:%d, right add:%d\n",countL,countR);	
}

//adding version,mzt is assigned inside
void RawData::GetXICData_dMET_NoLim(vector<CalData> &XICList,double mze,int ch,int MS2Scan,double mzt)
{	
	size_t ScanTotal=ScanMap.size();
	int MS1Scan=seekScan(MS2Scan,0,ScanTotal);
	int tempIDX=MS1Scan;
	CalData ct;
	double RTOld=MSList[MS1Scan].RT;
	ct.mzt=mzt;
	int interput=0;

	if(IS_DYA_XIC) DFT.ReInitial();

	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO_dMET(ct,ch,mzt,ppmLV,ppmHV);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;		
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{			
			XICList.push_back(ct);		
			interput=0;
			RTOld=MSList[tempIDX].RT;	
			if(IS_DYA_XIC) DFT.Add(-1*ct.RT,ct.GetmaxInt());
		}	
		
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT<MINRT+MSList[MS1Scan].RT) break;
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;
	interput=0;	
	RTOld=MSList[MS1Scan].RT;

	if(IS_DYA_XIC) DFT.ReInitial();

	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO_dMET(ct,ch,mzt,ppmLV,ppmHV);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;	
			else if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			XICList.push_back(ct);
			interput=0;
			RTOld=MSList[tempIDX].RT;	
			if(IS_DYA_XIC) DFT.Add(ct.RT,ct.GetmaxInt());
		}
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT>MAXRT+MSList[MS1Scan].RT) break;		
		tempIDX++;		
	}	
	//printf("Left add:%d, right add:%d\n",countL,countR);	
}

//adding version,mzt is assigned inside
void RawData::GetData_dMET(vector<CalData> &SList,double mze,int ch,int MS2Scan,double mzt)
{	
	size_t ScanTotal=ScanMap.size();
	int MS1Scan=seekScan(MS2Scan,0,ScanTotal);	
	CalData ct;
	double RTOld=MSList[MS1Scan].RT;
	ct.mzt=mzt;
	bool IsFind=MSList[MS1Scan].FindISO_dMET(ct,ch,mzt,ppmLV,ppmHV);
	if(IsFind) SList.push_back(ct);		
}

//assign mzt outside this function
void RawData::GetXICData(vector<CalData> &XICList,double pmz,int ch,int MS2Scan)
{
	XICList.clear();
	size_t ScanTotal=ScanMap.size();
	int MS1Scan=seekScan(MS2Scan,0,ScanTotal);
	int tempIDX=MS1Scan;
	CalData ct;
	double RTOld=MSList[MS1Scan].RT;
	int interput=0;

	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO(ct,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{
			XICList.push_back(ct);
			interput=0;
			RTOld=MSList[tempIDX].RT;
		}		
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;
	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO(ct,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			XICList.push_back(ct);
			interput=0;
			RTOld=MSList[tempIDX].RT;
		}
		
		tempIDX++;		
	}
}

//assign mzt outside this function
void RawData::GetData(vector<CalData> &SList,double pmz,int ch,int MS2Scan)
{
	SList.clear();
	size_t ScanTotal=ScanMap.size();
	int MS1Scan=seekScan(MS2Scan,0,ScanTotal);	
	CalData ct;
	double RTOld=MSList[MS1Scan].RT;
	bool IsFind=MSList[MS1Scan].FindISO(ct,ch,pmz);
	if(IsFind) SList.push_back(ct);
}

void RawData::RandomSelect(vector<CalData> &XICList)
{
	size_t i, total=XICList.size();
	if(total<=MAX_SVM_TRAIN) return;
	double Prob=MAX_SVM_TRAIN*1.0/total;
	vector<CalData> tmpdata;
	srand((unsigned)time(NULL));
	for(i=0;i<total;i++)
	{
		int r=rand();
		if(r>32767*Prob) continue;
		tmpdata.push_back(XICList[i]);		
	}
	total=tmpdata.size();
	XICList.resize(total);
	for(i=0;i<total;i++) XICList[i]=tmpdata[i];
}

bool RawData::ModelBuilding(vector<CalData> &XICList)
{
	if(CalModel!=NULL)
	{
		delete CalModel;
		CalModel=NULL;
	}
	CalModel=new SVMModel;
	return CalModel->Model_Train(XICList,SelParIdx);
	//double RS=CalModel->ParmSelection(XICList,SelParIdx,3);
	//if(RS<0.1) return false;
	//return true;
}

bool RawData::EstimateCalMET(vector<Pre_calData> &IDCalData)
{
	vector<CalData> XICList;	
	int pnum=IDCalData.size();	
	if(pnum<=0) return false;
	modeTotal=pnum;
	double *ppmdata;
	ppmdata=new double[pnum];	
	int k=0;
	for(int i=0;i<pnum;i++)
	{
		GetXICData(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
		if(XICList.size()==0) continue;		
		ppmdata[k]=(GetCalPmz(XICList)-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;	
		k++;
		XICList.clear();
	}
	mtB.Convert(ppmdata,k);	
	OutlierRem ot;
	if(ot.ORemMixModel(ppmdata,k))
	{
		mtB.mean=ot.m_new;
		mtB.std=ot.std_new;
	}
	else
	{
		mtB.mean=gsl_stats_mean(ppmdata,1,k);
		mtB.std=gsl_stats_sd_m(ppmdata,1,k,mtB.mean);
	}
	delete []ppmdata;	
	return true;
}

double RawData::GetMaxInt(vector<CalData> &XICList)
{
	size_t i,n=XICList.size();
	double maxInt=0;
	for(i=0;i<n;i++)
	{
		double tmpf=XICList[i].GetmaxIntR();
		if(maxInt<tmpf) maxInt=tmpf;
	}
	return maxInt;
}

bool RawData::EstimateCalMET_Lin(vector<Pre_calData> &IDCalData)
{
	vector<CalData> XICList;	
	int pnum=IDCalData.size();	
	if(pnum<=0) return false;
	modeTotal=pnum;
	double *ppmdata,*IntData;
	ppmdata=new double[pnum];
	IntData=new double[pnum];
	int k=0;
	for(int i=0;i<pnum;i++)
	{
		GetXICData_dMET_NoLim(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
		if(XICList.size()==0) continue;		
		ppmdata[k]=(GetCalPmz_Lin(XICList)-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;		
		IntData[k]=XICList.size();
		k++;
		XICList.clear();
	}
	mtB.Convert(ppmdata,k);	
	OutlierRem ot;
	if(ot.ORemMixModel(ppmdata,k))
	{
		mtB.mean=ot.m_new;
		mtB.std=ot.std_new;
	}
	else
	{
		mtB.mean=gsl_stats_mean(ppmdata,1,k);
		mtB.std=gsl_stats_sd_m(ppmdata,1,k,mtB.mean);
	}
	IsMETModel=false;
	printf("begin to modeling the intensity relative model:\n");
	if(IntR_ModelBuildingNew(ppmdata,IntData,k,METModel))IsMETModel=true;
	printf("end modeling\n");
	//for debug and test
	if(OIntR_Data)OutPutModelData(IntRModelData,ppmdata,IntData,k);
	//
	delete []IntData;
	delete []ppmdata;	
	return true;
}

bool RawData::ModelBuilding_Non(vector<Pre_calData> &IDCalData)
{
	GetPreMET(IDCalData);//calculation the original 
	vector<CalData> XICList;	
	int pnum=IDCalData.size();	
	if(pnum<=0) return false;
	modeTotal=pnum;
	double *ppmdata,*IntData;
	ppmdata=new double[pnum];
	IntData=new double[pnum];
	int k=0;
	for(int i=0;i<pnum;i++)
	{
		GetXICData_dMET_NoLim(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
		if(XICList.size()==0) continue;		
		ppmdata[k]=(GetCalPmz_Non(XICList)-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;	
		IntData[k]=XICList.size();
		k++;
		XICList.clear();
	}
	mtB.Convert(ppmdata,k);	
	OutlierRem ot;
	if(ot.ORemMixModel(ppmdata,k))
	{
		mtB.mean=ot.m_new;
		mtB.std=ot.std_new;
	}
	else
	{
		mtB.mean=gsl_stats_mean(ppmdata,1,k);
		mtB.std=gsl_stats_sd_m(ppmdata,1,k,mtB.mean);
	}
	IsMETModel=false;
	printf("begin to modeling the intensity relative model:\n");
	if(IntR_ModelBuildingNew(ppmdata,IntData,k,METModel))IsMETModel=true;
	printf("end modeling\n");
	//for debug and test
	if(OIntR_Data)OutPutModelData(IntRModelData,ppmdata,IntData,k);
	//
	delete []IntData;
	delete []ppmdata;	
	return true;	
}

bool RawData::EstimateCalMET_dMET(vector<Pre_calData> &IDCalData)
{
	vector<CalData> XICList;	
	int pnum=IDCalData.size();	
	if(pnum<=0) return false;
	modeTotal=pnum;
	double *ppmdata,*IntData;
	ppmdata=new double[pnum];
	IntData=new double[pnum];
	int k=0;
	for(int i=0;i<pnum;i++)
	{
		GetXICData_dMET_NoLim(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
		if(XICList.size()==0) continue;		
		ppmdata[k]=(GetCalPmz(XICList)-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;	
		IntData[k]=XICList.size();
		k++;
		XICList.clear();
	}
	mtB.Convert(ppmdata,k);	
	OutlierRem ot;
	if(ot.ORemMixModel(ppmdata,k))
	{
		mtB.mean=ot.m_new;
		mtB.std=ot.std_new;
	}
	else
	{
		mtB.mean=gsl_stats_mean(ppmdata,1,k);
		mtB.std=gsl_stats_sd_m(ppmdata,1,k,mtB.mean);
	}
	IsMETModel=false;
	printf("begin to modeling the intensity relative model:\n");
	if(IntR_ModelBuildingNew(ppmdata,IntData,k,METModel))IsMETModel=true;
	printf("end modeling\n");
	//for debug and test
	if(OIntR_Data)OutPutModelData(IntRModelData,ppmdata,IntData,k);
	//
	delete []IntData;
	delete []ppmdata;	
	return true;
}

bool RawData::EstimateCalMET_S(vector<Pre_calData> &IDCalData)
{
	vector<CalData> SList;	
	int pnum=IDCalData.size();	
	if(pnum<=0) return false;
	modeTotal=pnum;
	double *ppmdata,*IntData;
	ppmdata=new double[pnum];
	IntData=new double[pnum];
	int k=0;
	for(int i=0;i<pnum;i++)
	{
		GetData_dMET(SList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
		if(SList.size()==0) continue;		
		ppmdata[k]=(GetCalPmz(SList)-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;	
		IntData[k]=SList.size();
		k++;
		SList.clear();
	}
	mtB.Convert(ppmdata,k);	
	OutlierRem ot;
	if(ot.ORemMixModel(ppmdata,k))
	{
		mtB.mean=ot.m_new;
		mtB.std=ot.std_new;
	}
	else
	{
		mtB.mean=gsl_stats_mean(ppmdata,1,k);
		mtB.std=gsl_stats_sd_m(ppmdata,1,k,mtB.mean);
	}
	IsMETModel=false;
	printf("begin to modeling the intensity relative model:\n");
	if(IntR_ModelBuildingNew(ppmdata,IntData,k,METModel))IsMETModel=true;
	printf("end modeling\n");
	//for debug and test
	if(OIntR_Data)OutPutModelData(IntRModelData,ppmdata,IntData,k);
	//
	delete []IntData;
	delete []ppmdata;	
	return true;
}

bool RawData::EstimateCalMET_LinS(vector<Pre_calData> &IDCalData)
{
	vector<CalData> SList;	
	int pnum=IDCalData.size();	
	if(pnum<=0) return false;
	modeTotal=pnum;
	double *ppmdata,*IntData;
	ppmdata=new double[pnum];
	IntData=new double[pnum];
	int k=0;
	for(int i=0;i<pnum;i++)
	{
		GetData_dMET(SList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
		if(SList.size()==0) continue;		
		ppmdata[k]=(GetCalPmz_Lin(SList)-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;	
		IntData[k]=SList.size();
		k++;
		SList.clear();
	}
	mtB.Convert(ppmdata,k);	
	OutlierRem ot;
	if(ot.ORemMixModel(ppmdata,k))
	{
		mtB.mean=ot.m_new;
		mtB.std=ot.std_new;
	}
	else
	{
		mtB.mean=gsl_stats_mean(ppmdata,1,k);
		mtB.std=gsl_stats_sd_m(ppmdata,1,k,mtB.mean);
	}
	IsMETModel=false;
	printf("begin to modeling the intensity relative model:\n");
	if(IntR_ModelBuildingNew(ppmdata,IntData,k,METModel)) IsMETModel=true;
	printf("end modeling\n");
	//for debug and test
	if(OIntR_Data)OutPutModelData(IntRModelData,ppmdata,IntData,k);
	//
	delete []IntData;
	delete []ppmdata;	
	return true;
}

void RawData::OutPutCalData(vector<CalData> &CalData,FILE *fp)
{
	size_t sg=CalData.size();
	CString tmpStr;
	for(size_t i=0;i<sg;i++)
	{
		CalData[i].OutPut2Str(tmpStr,SelParIdx);
		fprintf(fp,"%s",tmpStr.GetBuffer());
		tmpStr.ReleaseBuffer();
	}
}

bool RawData::OutPutCalData(vector<Pre_calData> &IDCalData,char *fname)
{
	FILE *fp;
	fp=fopen(fname,"w");
	if(fp==NULL) return false;
	vector<CalData> XICList;	
	int pnum=IDCalData.size();	
	modeTotal=pnum;
	for(int i=0;i<pnum;i++)
	{
		GetXICData(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
	}
	OutPutCalData(XICList,fp);	
	return true;
}

void RawData::GetPreMET(vector<Pre_calData> &IDCalData)
{		
	size_t pnum=IDCalData.size();	
	if(pnum<=0)
	{
		ppmLV=-MET;
		ppmHV=MET;
		return;
	}
	vector<double> ppme;	
	for(size_t i=0;i<pnum;i++)
	{
		double ppmt=(IDCalData[i].mze-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;
		if(ppmt>PPMLIM_MIN&&ppmt<PPMLIM_MAX)ppme.push_back(ppmt);//a  safe range check
	}

	mtA.Convert(ppme);	
	OutlierRem ot;
	if(ot.ORemMixModel(ppme))
	{
		if(ot.std_new>1e-4)
		{
			ppmLV=ot.m_new-3*ot.std_new;
			ppmHV=ot.m_new+3*ot.std_new;
			mtA.mean=ot.m_new;
			mtA.std=ot.std_new;
			MET=3*ot.std_new*SQRTTWO;
		}
		else
		{
			ppmLV=-MET;
			ppmHV=MET;
		}
	}
	else
	{
		ppmLV=-MET;
		ppmHV=MET;
		mygsl_stats_stdm(ppme,&(mtA.mean),&(mtA.std));
	}	
}

bool RawData::ModelBuilding_dMET(vector<Pre_calData> &IDCalData)
{
	if(XICCALIBRATE) return ModelBuilding_X(IDCalData);
	else return ModelBuilding_S(IDCalData); 	
}


bool RawData::ModelBuilding_Lin(vector<Pre_calData> &IDCalData)
{
	if(XICCALIBRATE) return ModelBuilding_LinX(IDCalData);
	else return ModelBuilding_LinS(IDCalData); 	
}

bool RawData::ModelBuilding_X(vector<Pre_calData> &IDCalData)
{
	GetPreMET(IDCalData);//calculation the original 
	vector<CalData> XICList;	
	size_t pnum=IDCalData.size();	
	modeTotal=pnum;
	printf("Begin to build the svm model:\n" );
	clock_t start,finish; 
	start=clock();
	printf("Step 1: prepare the data.\n" );
	for(size_t i=0;i<pnum;i++)
	{
		double ppmAct=(IDCalData[i].mze-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;
		if(ppmAct>ppmHV||ppmAct<ppmLV) continue;
		GetXICData_dMET(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);		
	}	
	printf("End of step 1.\n");	
	printf("Step 2: select the case to train the svm model.\n");
	RandomSelect(XICList);		
	printf("End of step 2\n");
	printf("Step 3: select the feature to train the svm model.\n");
	featureSel(XICList);
	printf("End of step 3.\n");
	finish=clock();
	printf("The prepare work consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);	
	start=finish;
	printf("Step 4: begin to train the svm model.\n");
	if(!ModelBuilding(XICList))
	{
		//printf("failure to build the calibration model.\n");
		Error_code|=ERROR_MDFAILUR;		
		return false;
	}
	printf("End of step 4.\n");
	finish=clock();
	printf("Training consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	printf("Estimate the svm residue mean and sigma.\n");
	start=finish;
	EstimateCalMET_dMET(IDCalData);	
	finish=clock();
	printf("End of this step: %lf.\n",(double)(finish-start)/CLOCKS_PER_SEC);	
	printf("The suggested MET is %lf.\n",mtB.std*3);	
	if(IsOutput) 
	{
		printf("Output Model data.\n");
		OutPutCalData(XICList,fp_cal_data);
		printf("End of Output Model data.\n");
	}
	return true;
}


bool RawData::ModelBuilding_S(vector<Pre_calData> &IDCalData)
{
	GetPreMET(IDCalData);//calculation the original 
	vector<CalData> SList;	
	size_t pnum=IDCalData.size();	
	modeTotal=pnum;
	printf("Begin to build the svm model:\n" );
	clock_t start,finish; 
	start=clock();
	printf("Step 1: prepare the data.\n" );
	for(size_t i=0;i<pnum;i++)
	{
		double ppmAct=(IDCalData[i].mze-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;
		if(ppmAct>ppmHV||ppmAct<ppmLV) continue;
		GetData_dMET(SList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);		
	}
	printf("End of step 1.\n");	
	printf("Step 2: select the case to train the svm model.\n");
	RandomSelect(SList);	
	printf("End of step 2\n");
	printf("Step 3: select the feature to train the svm model.\n");
	featureSel(SList);
	printf("End of step 3.\n");
	finish=clock();
	printf("The prepare work consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);	
	start=finish;
	printf("Step 4: begin to train the svm model.\n");
	if(!ModelBuilding(SList))
	{
		//printf("failure to build the calibration model.\n");
		Error_code|=ERROR_MDFAILUR;		
		return false;
	}
	printf("End of step 4.\n");
	finish=clock();
	printf("Training consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	printf("Estimate the svm residue mean and sigma.\n");
	start=finish;
	EstimateCalMET_S(IDCalData);
	finish=clock();
	printf("End of this step: %lf.\n",(double)(finish-start)/CLOCKS_PER_SEC);
	printf("The suggested MET is %lf.\n",mtB.std*3);
	if(IsOutput)
	{
		printf("Output Model data.\n");
		OutPutCalData(SList,fp_cal_data);
		printf("End of Output Model data.\n");		
	}
	return true;
}

bool RawData::ModelBuilding_LinX(vector<Pre_calData> &IDCalData)
{
	GetPreMET(IDCalData);//calculation the original 
	vector<CalData> XICList;	
	size_t pnum=IDCalData.size();	
	modeTotal=pnum;
	printf("Begin to build the linear model:\n" );
	clock_t start,finish; 
	start=clock();
	printf("Step 1: prepare the data.\n" );
	for(size_t i=0;i<pnum;i++)
	{
		double ppmAct=(IDCalData[i].mze-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;
		if(ppmAct>ppmHV||ppmAct<ppmLV) continue;
		GetXICData_dMET(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);		
	}	
	printf("End of step 1.\n");	
	printf("Step 2: select the case to train the linear model.\n");
	RandomSelect(XICList);		
	printf("End of step 2\n");
	printf("Step 3: select the feature to train the linear model.\n");
	featureSel(XICList);
	printf("End of step 3.\n");
	finish=clock();
	printf("The prepare work consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);	
	start=finish;
	printf("Step 4: begin to regress the linear model.\n");
	if(!LModel.ModelTrain(XICList,SelParIdx))
	{
		//printf("failure to build the calibration model.\n");
		Error_code|=ERROR_MDFAILUR;		
		return false;
	}
	printf("End of step 4.\n");
	finish=clock();
	printf("Training consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	printf("Estimate the svm residue mean and sigma.\n");
	start=finish;
	EstimateCalMET_Lin(IDCalData);	
	finish=clock();
	printf("End of this step: %lf.\n",(double)(finish-start)/CLOCKS_PER_SEC);	
	printf("The suggested MET is %lf.\n",mtB.std*3);	
	if(IsOutput) 
	{
		printf("Output Model data.\n");
		OutPutCalData(XICList,fp_cal_data);
		printf("End of Output Model data.\n");
	}
	return true;
}


bool RawData::ModelBuilding_LinS(vector<Pre_calData> &IDCalData)
{
	GetPreMET(IDCalData);//calculation the original 
	vector<CalData> SList;	
	size_t pnum=IDCalData.size();	
	modeTotal=pnum;
	printf("Begin to build the linear model:\n" );
	clock_t start,finish; 
	start=clock();
	printf("Step 1: prepare the data.\n" );
	for(size_t i=0;i<pnum;i++)
	{
		double ppmAct=(IDCalData[i].mze-IDCalData[i].mzt)*1e6/IDCalData[i].mzt;
		if(ppmAct>ppmHV||ppmAct<ppmLV) continue;
		GetData_dMET(SList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);		
	}
	printf("End of step 1.\n");	
	printf("Step 2: select the case to train the linear model.\n");
	RandomSelect(SList);	
	printf("End of step 2\n");
	printf("Step 3: select the feature to train the linear model.\n");
	featureSel(SList);
	printf("End of step 3.\n");
	finish=clock();
	printf("The prepare work consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);	
	start=finish;
	printf("Step 4: begin to train the linear model.\n");
	if(!LModel.ModelTrain(SList,SelParIdx))
	{
		//printf("failure to build the calibration model.\n");
		Error_code|=ERROR_MDFAILUR;		
		return false;
	}
	printf("End of step 4.\n");
	finish=clock();
	printf("Training consume time: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	printf("Estimate the linear residue mean and sigma.\n");
	start=finish;
	EstimateCalMET_LinS(IDCalData);
	finish=clock();
	printf("End of this step: %lf.\n",(double)(finish-start)/CLOCKS_PER_SEC);
	printf("The suggested MET is %lf.\n",mtB.std*3);
	if(IsOutput)
	{
		printf("Output Model data.\n");
		OutPutCalData(SList,fp_cal_data);
		printf("End of Output Model data.\n");		
	}
	return true;
}


void RawData::GetPreHis(vector<CalData> &XICList)
{
	size_t i,n=XICList.size();
	if(n<=0) return;
	double *ppm;
	ppm=new double[n];	
	for(i=0;i<n;i++)ppm[i]=(XICList[i].Isomz[0]-XICList[i].mzt)*1e6/XICList[i].mzt;	
	mtA.Convert(ppm,n);	
	OutlierRem ot;
	ot.ORemMixModel(ppm,n);
	mtA.mean=ot.m_new;
	mtA.std=ot.std_new;
	delete []ppm;	
}

void RawData::GetCalHis(vector<CalData> &XICList)
{
	if(CalModel==NULL) return;
	size_t i,n=XICList.size();
	if(n<=0) return;
	double *ppm;
	ppm=new double[n];	
	for(i=0;i<n;i++)
	{
		double mzp=CalModel->Predict(XICList[i],SelParIdx);	
		ppm[i]=(mzp-XICList[i].mzt)*1e6/XICList[i].mzt;	
	}
	mtB.Convert(ppm,n);	
	OutlierRem ot;
	ot.ORemMixModel(ppm,n);
	mtB.mean=ot.m_new;
	mtB.std=ot.std_new;
	delete []ppm;	
}

void RawData::GetCalHis(vector<CalData> &XICList,double bT[2])
{
	if(CalModel==NULL) return;
	size_t i,n=XICList.size();
	if(n<=0) return;
	double *ppm;
	ppm=new double[n];	
	for(i=0;i<n;i++)
	{
		double mzp=CalModel->Predict(XICList[i],SelParIdx);	
		ppm[i]=(mzp-XICList[i].mzt)*1e6/XICList[i].mzt;	
	}
	OutlierRem ot;
	if(ot.oRemoveGTest(ppm,n))
	{
		bT[0]=ot.m_new;
		bT[1]=ot.std_new;
	}
	else
	{
		bT[0]=gsl_stats_mean(ppm,1,n);
		bT[1]=gsl_stats_sd_m(ppm,1,n,bT[0]);
	}
	delete []ppm;	
}

double RawData::GetCalPmz(vector<CalData> &XICList)
{
	if(CalModel==NULL) return 0;
	size_t i,n=XICList.size();
	if(n<=0) return 0;
	vector<double> pmzCal;
	vector<double> weight;
	for(i=0;i<n;i++)CalModel->Predict(XICList[i],pmzCal,weight,SelParIdx);
	OutlierRem ot;
	double bRetrun=0;
	if(ot.oRemoveGTestW(pmzCal,weight))bRetrun=ot.m_new;
	else bRetrun=ot.MeanW(pmzCal,weight);	
	return bRetrun;
}

bool RawData::GetCalPmz(vector<CalData> &XICList,double bRT[2])
{
	if(CalModel==NULL) return false;
	size_t i,n=XICList.size();
	if(n<=0) return false;
	vector<double>pmzCal;
	vector<double> weight;
	for(i=0;i<n;i++)CalModel->Predict(XICList[i],pmzCal,weight,SelParIdx);	
	if(pmzCal.size()<=0) return false;
	OutlierRem ot;
	double bRetrun=0;
	if(ot.oRemoveGTestW(pmzCal,weight))
	{
		bRT[0]=ot.m_new;
		bRT[1]=ot.std_new;
	}
	else 
	{
		bRT[0]=ot.MeanW(pmzCal,weight);
		bRT[1]=ot.StdW(bRT[0],pmzCal,weight);
	}	
	return true;
}

double RawData::GetCalPmz_Lin(vector<CalData> &XICList)
{
	if(LModel.stats[2]<0.1) return 0;
	size_t i,n=XICList.size();
	if(n<=0) return 0;
	vector<double> pmzCal;
	vector<double> weight;
	for(i=0;i<n;i++)LModel.Predict(XICList[i],pmzCal,weight,SelParIdx);


	OutlierRem ot;
	double bRetrun=0;
	if(ot.oRemoveGTestW(pmzCal,weight))bRetrun=ot.m_new;
	else bRetrun=ot.MeanW(pmzCal,weight);	
	return bRetrun;
}

double RawData::GetCalPmz_Non(vector<CalData> &XICList)
{
	size_t i,n=XICList.size();
	if(n<=0) return 0;
	vector<double> pmzCal;
	vector<double> weight;
	for(i=0;i<n;i++)
	{
		for(int k=0;k<XICList[i].IsoNum;k++)
		{
			pmzCal.push_back(XICList[i].Isomz[i]-ISO_DIFF[i]/XICList[i].ch);
			weight.push_back(sqrt(XICList[i].IsotopicE[i]));
		}
	}
	OutlierRem ot;
	double bRetrun=0;
	if(ot.oRemoveGTestW(pmzCal,weight))bRetrun=ot.m_new;
	else bRetrun=ot.MeanW(pmzCal,weight);	
	return bRetrun;
}

bool RawData::GetCalPmz_Lin(vector<CalData> &XICList,double bRT[2])
{
	if(LModel.stats[2]<0.1) return 0;
	size_t i,n=XICList.size();
	if(n<=0) return false;
	vector<double>pmzCal;
	vector<double> weight;
	for(i=0;i<n;i++)LModel.Predict(XICList[i],pmzCal,weight,SelParIdx);
	
	if(pmzCal.size()<=0) return false;
	OutlierRem ot;
	double bRetrun=0;
	if(ot.oRemoveGTestW(pmzCal,weight))
	{
		bRT[0]=ot.m_new;
		bRT[1]=ot.std_new;
	}
	else 
	{
		bRT[0]=ot.MeanW(pmzCal,weight);
		bRT[1]=ot.StdW(bRT[0],pmzCal,weight);
	}	
	return true;
}

bool RawData::GetCalPmz_Non(vector<CalData> &XICList,double bRT[2])
{	
	size_t i,n=XICList.size();
	if(n<=0) return false;
	vector<double>pmzCal;
	vector<double> weight;
	for(i=0;i<n;i++)
	{
		for(int k=0;k<XICList[i].IsoNum;k++)
		{
			pmzCal.push_back(XICList[i].Isomz[i]-ISO_DIFF[i]/XICList[i].ch);
			weight.push_back(sqrt(XICList[i].IsotopicE[i]));
		}

	}	
	if(pmzCal.size()<=0) return false;
	OutlierRem ot;
	double bRetrun=0;
	if(ot.oRemoveGTestW(pmzCal,weight))
	{
		bRT[0]=ot.m_new;
		bRT[1]=ot.std_new;
	}
	else 
	{
		bRT[0]=ot.MeanW(pmzCal,weight);
		bRT[1]=ot.StdW(bRT[0],pmzCal,weight);
	}	
	return true;
}


int RawData::ProcessOneXICBack(CalData &ct,int MS1Scan,double bReturn[2])
{
	if(CalModel==NULL) return 0;
	int tempIDX=MS1Scan-1;
	vector<double>ppmD;
	double RTOld=MSList[MS1Scan].RT;	
	CalModel->Predict(ct,ppmD,SelParIdx);	
	double pmz=ct.Isomz[0];
	int ch=ct.ch;
	int interput=0;
	CalData tmpCT;
	int ScanTotal=MSList.size();
	int bCount=0;
	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{
			CalModel->Predict(tmpCT,ppmD,SelParIdx);	
			bCount++;		
			interput=0;
			RTOld=MSList[tempIDX].RT;
		}
		
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;
	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			CalModel->Predict(tmpCT,ppmD,SelParIdx);
			bCount++;			
			interput=0;
			RTOld=MSList[tempIDX].RT;
		}
		
		tempIDX++;		
	}
	if(ppmD.size()<=0) return 0;	
	OutlierRem ot;
	if(ot.ORemMixModel(ppmD))
	{
		bReturn[0]=ot.m_new*1e6/(mtB.mean+1e6);
		bReturn[1]=ot.std_new;
	}
	else
	{
		bReturn[0]=ot.Mean(ppmD)*1e6/(1e6+mtB.mean);
		bReturn[1]=ot.Std(bReturn[0],ppmD);
	}	
	return bCount;
}

//process the known charge case
int RawData::ProcessOneXIC(CalData &ct,int MS1Scan,double bReturn[2])
{
	if(CalModel==NULL) return 0;
	int tempIDX=MS1Scan-1;
	double RTOld=MSList[MS1Scan].RT;
	vector<CalData> PCTList;
	PCTList.push_back(ct);
	double pmz=ct.Isomz[0];
	int ch=ct.ch;
	int interput=0;
	CalData tmpCT;
	int ScanTotal=MSList.size();
	int bCount=0;

	if(IS_DYA_XIC) DFT.ReInitial();

	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{
			PCTList.push_back(tmpCT);
			interput=0;
			RTOld=MSList[tempIDX].RT;
			if(IS_DYA_XIC) DFT.Add(-1*tmpCT.RT,tmpCT.GetmaxInt());
		}	
		//added on 2012.2.11
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT<MINRT+MSList[MS1Scan].RT) break;
		//end
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;
	if(IS_DYA_XIC) DFT.ReInitial();
	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			PCTList.push_back(tmpCT);
			interput=0;
			RTOld=MSList[tempIDX].RT;
			if(IS_DYA_XIC) DFT.Add(tmpCT.RT,tmpCT.GetmaxInt());
		}
		//added on 2012.2.11
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT>MAXRT+MSList[MS1Scan].RT) break;
		//end		
		tempIDX++;		
	}
	if(!GetCalPmz(PCTList,bReturn)) return 0;
	bReturn[0]=bReturn[0]*1e6/(mtB.mean+1e6);
	return PCTList.size();
}


//process the known charge case
int RawData::ProcessOneXIC_Lin(CalData &ct,int MS1Scan,double bReturn[2])
{
	if(LModel.stats[2]<0.1) return 0;
	int tempIDX=MS1Scan-1;
	double RTOld=MSList[MS1Scan].RT;
	vector<CalData> PCTList;
	PCTList.push_back(ct);
	double pmz=ct.Isomz[0];
	int ch=ct.ch;
	int interput=0;
	CalData tmpCT;
	int ScanTotal=MSList.size();
	int bCount=0;

	if(IS_DYA_XIC) DFT.ReInitial();

	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{
			PCTList.push_back(tmpCT);
			interput=0;
			RTOld=MSList[tempIDX].RT;
			if(IS_DYA_XIC) DFT.Add(-1*tmpCT.RT,tmpCT.GetmaxInt());
		}	
		//added on 2012.2.11
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT<MINRT+MSList[MS1Scan].RT) break;
		//end
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;
	if(IS_DYA_XIC) DFT.ReInitial();
	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			PCTList.push_back(tmpCT);
			interput=0;
			RTOld=MSList[tempIDX].RT;
			if(IS_DYA_XIC) DFT.Add(tmpCT.RT,tmpCT.GetmaxInt());
		}
		//added on 2012.2.11
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT>MAXRT+MSList[MS1Scan].RT) break;
		//end		
		tempIDX++;		
	}
	if(!GetCalPmz_Lin(PCTList,bReturn)) return 0;
	bReturn[0]=bReturn[0]*1e6/(mtB.mean+1e6);
	return PCTList.size();
}

//process the known charge case
int RawData::ProcessOneXIC_Non(CalData &ct,int MS1Scan,double bReturn[2])
{	
	int tempIDX=MS1Scan-1;
	double RTOld=MSList[MS1Scan].RT;
	vector<CalData> PCTList;
	PCTList.push_back(ct);
	double pmz=ct.Isomz[0];
	int ch=ct.ch;
	int interput=0;
	CalData tmpCT;
	int ScanTotal=MSList.size();
	int bCount=0;

	if(IS_DYA_XIC) DFT.ReInitial();

	while(tempIDX>=0)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT+MAX_INT_RT<RTOld) break;
		}
		else 
		{
			PCTList.push_back(tmpCT);
			interput=0;
			RTOld=MSList[tempIDX].RT;
			if(IS_DYA_XIC) DFT.Add(-1*tmpCT.RT,tmpCT.GetmaxInt());
		}	
		//added on 2012.2.11
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT<MINRT+MSList[MS1Scan].RT) break;
		//end
		tempIDX--;		
	}

	tempIDX=MS1Scan+1;
	if(IS_DYA_XIC) DFT.ReInitial();
	while(tempIDX<ScanTotal)
	{
		bool IsFind=MSList[tempIDX].FindISO(tmpCT,ch,pmz);
		if(!IsFind) 
		{
			interput++;
			if(interput>=INT_TIME_MAX) break;
			if(MSList[tempIDX].RT>MAX_INT_RT+RTOld) break;
		}
		else 
		{
			PCTList.push_back(tmpCT);
			interput=0;
			RTOld=MSList[tempIDX].RT;
			if(IS_DYA_XIC) DFT.Add(tmpCT.RT,tmpCT.GetmaxInt());
		}
		//added on 2012.2.11
		if(IS_DYA_XIC)
		{
			if(DFT.IsMinmal()) break;
		}
		else if(MSList[tempIDX].RT>MAXRT+MSList[MS1Scan].RT) break;
		//end		
		tempIDX++;		
	}
	if(!GetCalPmz_Non(PCTList,bReturn)) return 0;
	bReturn[0]=bReturn[0]*1e6/(mtB.mean+1e6);
	return PCTList.size();
}


bool RawData::XICCalibrate(size_t MS2,vector<CalReturn>& bReturn)
{
	if(PARIONRED) return XICCalibrateC(MS2,bReturn);
	else return XICCalibrateS(MS2,bReturn);
}

bool RawData::XICCalibrate_Lin(size_t MS2,vector<CalReturn>& bReturn)
{
	if(PARIONRED) return XICCalibrateC_Lin(MS2,bReturn);
	else return XICCalibrateS_Lin(MS2,bReturn);
}

bool RawData::XICCalibrate_Non(size_t MS2,vector<CalReturn>& bReturn)
{
	if(PARIONRED) return XICCalibrateC_Non(MS2,bReturn);
	else return XICCalibrateS_Non(MS2,bReturn);
}

bool RawData::XICCalibrateC(size_t MS2,vector<CalReturn>& bReturn)
{	
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;	
	int MS1Scan=tmpIDX;	
	vector<CalData> PIso;
	//at first, we should find all the possible parent ions
	string stmp;
	char buf[256];
	stmp="Begin to Re-determinate the parent ions:\n";
	DInfBuf.push_back(stmp);

	clock_t start,finish; 
	start=clock();

#ifdef DEBUG_CASE
	//for raw file B06-11071.RAW in dir: G:\BNPModelForMascot\FTcontrolDataset\raw
	if(MS2Scans[MS2].scan==5129)
	{
		int case_i=0;//no meaning, just to stop the process here
		case_i++;
	}
#endif

	MSList[MS1Scan].ExtendFindISO(PIso,MS2Scans[MS2].pmz,MS2Scans[MS2].IW);

#ifdef DEBUG_CASE
	//for raw file B06-11071.RAW in dir: G:\BNPModelForMascot\FTcontrolDataset\raw
	if(MS2Scans[MS2].scan==5129)
	{
		int case_i=0;//no meaning, just to stop the process here
		case_i++;
	}
#endif

	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);
	start=finish;
	sprintf(buf,"Begin to search the XIC data and calibrate for %d parent ions.\n",PIso.size());
	stmp=buf;
	DInfBuf.push_back(stmp);

#ifdef CASE_TEST
	int pionnum=PIso.size();
	fprintf(ParentIons,"%d\t%d\t%d\t%lf\n",pionnum,MS2Scans[MS2].scan,MS2Scans[MS2].ch,MS2Scans[MS2].pmz);	
#endif

	size_t i, pPNum=PIso.size();
	if(pPNum<=0) return false;
	double bR[2]; 
	vector<CalData> tmpPISO;
	CalReturn crt;
	
	for(i=0;i<pPNum;i++)
	{		
		if(PIso[i].ch!=CH_UNKNOWN)
		{
			crt.spectrum_num=ProcessOneXIC(PIso[i],MS1Scan,bR);
			if(crt.spectrum_num>0)
			{
				crt.cal_mean=bR[0];
				crt.cal_std=bR[1];
				crt=PIso[i];
				bReturn.push_back(crt);
			}
		}
		else//尝试利用周围图谱确定电荷
		{
			CalData CT=PIso[i];
			vector<int> CHL;
			GetPossChFromPmz(PIso[i].Isomz[0],PIso[i].scannum,CHL);
			size_t sg=CHL.size();
			if(sg>0)
			{
				for(size_t k=0;k<sg;k++)
				{
					CT.ch=CHL[k];
					tmpPISO.push_back(CT);
				}
			}
			else //实在无法确定电荷，就以电荷不确定论处
			{
				crt.spectrum_num=ProcessOneXIC(PIso[i],MS1Scan,bR);
				if(crt.spectrum_num>0)
				{
					crt.cal_mean=bR[0];
					crt.cal_std=bR[1];
					crt=PIso[i];					
					bReturn.push_back(crt);
				}
			}
		}
	}

	pPNum=tmpPISO.size();//补充确定周围电荷的情况
	for(i=0;i<pPNum;i++)
	{
		crt.spectrum_num=ProcessOneXIC(tmpPISO[i],MS1Scan,bR);
		if(crt.spectrum_num>0)
		{
			crt.cal_mean=bR[0];
			crt.cal_std=bR[1];
			crt=tmpPISO[i];		
			bReturn.push_back(crt);
		}
	}
	finish=clock();
	sprintf(buf,"End to search the XIC data and calibrate: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);	
	return true;
}

bool RawData::XICCalibrateS(size_t MS2,vector<CalReturn>& bReturn)
{
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;	
	string stmp;
	char buf[256];
	stmp="Begin to perform the XIC calibration:\n";
	DInfBuf.push_back(stmp);
	clock_t start,finish; 
	start=clock();
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;
	int MS1Scan=tmpIDX;	
	CalData PIso;
	CalReturn ct;
	double bR[2];
	bool IsFind=false;
	//at first, we should find all the possible parent ions
	if(MSList[MS1Scan].FindISO(PIso,MS2Scans[MS2].ch,MS2Scans[MS2].pmz))
	{
		ct.spectrum_num=ProcessOneXIC(PIso,MS1Scan,bR);
		if(ct.spectrum_num>0)
		{
			ct.cal_mean=bR[0];
			ct.cal_std=bR[1];
			ct=PIso;
			bReturn.push_back(ct);
			IsFind=true;
		}
	}
	if(!IsFind)
	{
		ct.charge=MS2Scans[MS2].ch;
		ct.cal_mean=MS2Scans[MS2].pmz-(ppmHV+ppmLV)/2;
		ct.iso_goodness=0.0;	
		ct.iso_num=0;
		ct.rel_int=0;
		ct.abs_int=0;
		ct.cal_std=(ppmHV-ppmLV)/6;
		bReturn.push_back(ct);
	}
	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);
	return true;
}


bool RawData::XICCalibrateC_Lin(size_t MS2,vector<CalReturn>& bReturn)
{	
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;	
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;	
	int MS1Scan=tmpIDX;	
	vector<CalData> PIso;
	//at first, we should find all the possible parent ions
	string stmp;
	char buf[256];
	stmp="Begin to Re-determinate the parent ions:\n";
	DInfBuf.push_back(stmp);

	clock_t start,finish; 
	start=clock();

#ifdef DEBUG_CASE
	//for raw file B06-11071.RAW in dir: G:\BNPModelForMascot\FTcontrolDataset\raw
	if(MS2Scans[MS2].scan==5129)
	{
		int case_i=0;//no meaning, just to stop the process here
		case_i++;
	}
#endif

	MSList[MS1Scan].ExtendFindISO(PIso,MS2Scans[MS2].pmz,MS2Scans[MS2].IW);

#ifdef DEBUG_CASE
	//for raw file B06-11071.RAW in dir: G:\BNPModelForMascot\FTcontrolDataset\raw
	if(MS2Scans[MS2].scan==5129)
	{
		int case_i=0;//no meaning, just to stop the process here
		case_i++;
	}
#endif

	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);
	start=finish;
	sprintf(buf,"Begin to search the XIC data and calibrate for %d parent ions.\n",PIso.size());
	stmp=buf;
	DInfBuf.push_back(stmp);

#ifdef CASE_TEST
	int pionnum=PIso.size();
	fprintf(ParentIons,"%d\t%d\t%d\t%lf\n",pionnum,MS2Scans[MS2].scan,MS2Scans[MS2].ch,MS2Scans[MS2].pmz);	
#endif

	size_t i, pPNum=PIso.size();
	if(pPNum<=0) return false;
	double bR[2]; 
	vector<CalData> tmpPISO;
	CalReturn crt;
	
	for(i=0;i<pPNum;i++)
	{	
		if(PIso[i].ch!=CH_UNKNOWN)
		{
			crt.spectrum_num=ProcessOneXIC_Lin(PIso[i],MS1Scan,bR);
			if(crt.spectrum_num>0)
			{
				crt.cal_mean=bR[0];
				crt.cal_std=bR[1];	
				crt=PIso[i];
				bReturn.push_back(crt);
			}
		}
		else//尝试利用周围图谱确定电荷
		{
			CalData CT=PIso[i];
			vector<int> CHL;
			GetPossChFromPmz(PIso[i].Isomz[0],PIso[i].scannum,CHL);
			size_t sg=CHL.size();
			if(sg>0)
			{
				for(size_t k=0;k<sg;k++)
				{
					CT.ch=CHL[k];
					tmpPISO.push_back(CT);
				}
			}
			else //实在无法确定电荷，就以电荷不确定论处
			{
				crt.spectrum_num=ProcessOneXIC_Lin(PIso[i],MS1Scan,bR);
				if(crt.spectrum_num>0)
				{
					crt.cal_mean=bR[0];
					crt.cal_std=bR[1];
					crt=PIso[i];					
					bReturn.push_back(crt);
				}
			}
		}
	}

	pPNum=tmpPISO.size();//补充确定周围电荷的情况
	for(i=0;i<pPNum;i++)
	{
		crt.spectrum_num=ProcessOneXIC_Lin(tmpPISO[i],MS1Scan,bR);
		if(crt.spectrum_num>0)
		{
			crt.cal_mean=bR[0];
			crt.cal_std=bR[1];
			crt=tmpPISO[i];		
			bReturn.push_back(crt);
		}
	}
	finish=clock();
	sprintf(buf,"End to search the XIC data and calibrate: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);	
	return true;
}

bool RawData::XICCalibrateS_Lin(size_t MS2,vector<CalReturn>& bReturn)
{
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;	
	string stmp;
	char buf[256];
	stmp="Begin to perform the XIC calibration:\n";
	DInfBuf.push_back(stmp);
	clock_t start,finish; 
	start=clock();
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;
	int MS1Scan=tmpIDX;	
	CalData PIso;
	CalReturn ct;
	double bR[2];
	bool IsFind=false;
	//at first, we should find all the possible parent ions
	if(MSList[MS1Scan].FindISO(PIso,MS2Scans[MS2].ch,MS2Scans[MS2].pmz))
	{
		ct.spectrum_num=ProcessOneXIC(PIso,MS1Scan,bR);
		if(ct.spectrum_num>0)
		{
			ct.cal_mean=bR[0];
			ct.cal_std=bR[1];
			ct=PIso;
			bReturn.push_back(ct);
			IsFind=true;
		}
	}
	if(!IsFind)
	{
		ct.charge=MS2Scans[MS2].ch;
		ct.cal_mean=MS2Scans[MS2].pmz-(ppmHV+ppmLV)/2;
		ct.iso_goodness=0.0;	
		ct.iso_num=0;
		ct.rel_int=0;
		ct.abs_int=0;
		ct.cal_std=(ppmHV-ppmLV)/6;
		bReturn.push_back(ct);
	}

	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);
	return true;
}



bool RawData::XICCalibrateC_Non(size_t MS2,vector<CalReturn>& bReturn)
{	
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;	
	int MS1Scan=tmpIDX;	
	vector<CalData> PIso;
	//at first, we should find all the possible parent ions
	string stmp;
	char buf[256];
	stmp="Begin to Re-determinate the parent ions:\n";
	DInfBuf.push_back(stmp);

	clock_t start,finish; 
	start=clock();

#ifdef DEBUG_CASE
	//for raw file B06-11071.RAW in dir: G:\BNPModelForMascot\FTcontrolDataset\raw
	if(MS2Scans[MS2].scan==5129)
	{
		int case_i=0;//no meaning, just to stop the process here
		case_i++;
	}
#endif

	MSList[MS1Scan].ExtendFindISO(PIso,MS2Scans[MS2].pmz,MS2Scans[MS2].IW);

#ifdef DEBUG_CASE
	//for raw file B06-11071.RAW in dir: G:\BNPModelForMascot\FTcontrolDataset\raw
	if(MS2Scans[MS2].scan==5129)
	{
		int case_i=0;//no meaning, just to stop the process here
		case_i++;
	}
#endif

	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);
	start=finish;
	sprintf(buf,"Begin to search the XIC data and calibrate for %d parent ions.\n",PIso.size());
	stmp=buf;
	DInfBuf.push_back(stmp);

#ifdef CASE_TEST
	int pionnum=PIso.size();
	fprintf(ParentIons,"%d\t%d\t%d\t%lf\n",pionnum,MS2Scans[MS2].scan,MS2Scans[MS2].ch,MS2Scans[MS2].pmz);	
#endif

	size_t i, pPNum=PIso.size();
	if(pPNum<=0) return false;
	double bR[2]; 
	vector<CalData> tmpPISO;
	CalReturn crt;
	
	for(i=0;i<pPNum;i++)
	{	
		if(PIso[i].ch!=CH_UNKNOWN)
		{
			crt.spectrum_num=ProcessOneXIC_Non(PIso[i],MS1Scan,bR);
			if(crt.spectrum_num>0)
			{
				crt.cal_mean=bR[0];
				crt.cal_std=bR[1];
				crt=PIso[i];
				bReturn.push_back(crt);
			}
		}
		else//尝试利用周围图谱确定电荷
		{
			CalData CT=PIso[i];
			vector<int> CHL;
			GetPossChFromPmz(PIso[i].Isomz[0],PIso[i].scannum,CHL);
			size_t sg=CHL.size();
			if(sg>0)
			{
				for(size_t k=0;k<sg;k++)
				{
					CT.ch=CHL[k];
					tmpPISO.push_back(CT);
				}
			}
			else //实在无法确定电荷，就以电荷不确定论处
			{
				crt.spectrum_num=ProcessOneXIC_Non(PIso[i],MS1Scan,bR);
				if(crt.spectrum_num>0)
				{
					crt.cal_mean=bR[0];
					crt.cal_std=bR[1];
					crt=PIso[i];					
					bReturn.push_back(crt);
				}
			}
		}
	}

	pPNum=tmpPISO.size();//补充确定周围电荷的情况
	for(i=0;i<pPNum;i++)
	{
		crt.spectrum_num=ProcessOneXIC_Non(tmpPISO[i],MS1Scan,bR);
		if(crt.spectrum_num>0)
		{
			crt.cal_mean=bR[0];
			crt.cal_std=bR[1];
			crt=tmpPISO[i];		
			bReturn.push_back(crt);
		}
	}
	finish=clock();
	sprintf(buf,"End to search the XIC data and calibrate: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);	
	return true;
}

bool RawData::XICCalibrateS_Non(size_t MS2,vector<CalReturn>& bReturn)
{
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;	
	string stmp;
	char buf[256];
	stmp="Begin to perform the XIC calibration:\n";
	DInfBuf.push_back(stmp);
	clock_t start,finish; 
	start=clock();
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;
	int MS1Scan=tmpIDX;	
	CalData PIso;
	CalReturn ct;
	double bR[2];
	bool IsFind=false;
	//at first, we should find all the possible parent ions
	if(MSList[MS1Scan].FindISO(PIso,MS2Scans[MS2].ch,MS2Scans[MS2].pmz))
	{
		ct.spectrum_num=ProcessOneXIC(PIso,MS1Scan,bR);
		if(ct.spectrum_num>0)
		{
			ct.cal_mean=bR[0];
			ct.cal_std=bR[1];
			ct=PIso;
			bReturn.push_back(ct);
			IsFind=true;
		}
	}
	if(!IsFind)
	{
		ct.charge=MS2Scans[MS2].ch;
		ct.cal_mean=MS2Scans[MS2].pmz-(ppmHV+ppmLV)/2;
		ct.iso_goodness=0.0;	
		ct.iso_num=0;
		ct.rel_int=0;
		ct.abs_int=0;
		ct.cal_std=(ppmHV-ppmLV)/6;
		bReturn.push_back(ct);
	}
	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions: %lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);
	return true;
}


bool RawData::SCalibrate(size_t MS2,vector<CalReturn>& bReturn)
{	
	if(PARIONRED) return SCalibrateC(MS2,bReturn);
	else return SCalibrateS(MS2,bReturn);	
}
bool RawData::SCalibrate_Lin(size_t MS2,vector<CalReturn>& bReturn)
{	
	if(PARIONRED) return SCalibrateC_Lin(MS2,bReturn);
	else return SCalibrateS_Lin(MS2,bReturn);	
}

bool RawData::SCalibrateC(size_t MS2,vector<CalReturn>& bReturn)
{	
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;
	int MS1Scan=tmpIDX;	
	vector<CalData> PIso;
	//at first, we should find all the possible parent ions
	char buf[256];
	string stmp="Begin to Re-determinate the parent ions:\n";
	DInfBuf.push_back(stmp);

	clock_t start,finish; 
	start=clock();
	MSList[MS1Scan].ExtendFindISO(PIso,MS2Scans[MS2].pmz,MS2Scans[MS2].IW);
	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions:%lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);

	
#ifdef CASE_TEST
	int pionnum=PIso.size();
	fprintf(ParentIons,"%d\t%d\t%d\t%lf\n",pionnum,MS2Scans[MS2].scan,MS2Scans[MS2].ch,MS2Scans[MS2].pmz);	
#endif

	size_t pPNum=PIso.size();
	if(pPNum<=0) return false;
	vector<double>ppmD;
	CalReturn crt;
	for(size_t i=0;i<pPNum;i++)
	{		
		CalModel->Predict(PIso[i],ppmD,SelParIdx);
		OutlierRem ot;
		if(ot.ORemMixModel(ppmD))
		{
			crt.cal_mean=ot.m_new;
			crt.cal_std=ot.std_new;
		}
		else
		{
			crt.cal_mean=ot.Mean(ppmD);
			crt.cal_std=ot.Std(crt.cal_mean,ppmD);
		}
		crt=PIso[i];
		bReturn.push_back(crt);			
		ppmD.clear();
	}
	return true;
}

bool RawData::SCalibrateS(size_t MS2,vector<CalReturn>& bReturn)
{	
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;
	int MS1Scan=tmpIDX;	
	CalData PIso;
	
	CalReturn crt;
	//at first, we should find all the possible parent ions
	if(MSList[MS1Scan].FindISO(PIso,MS2Scans[MS2].ch,MS2Scans[MS2].pmz))
	{
		vector<double>ppmD;		
		CalModel->Predict(PIso,ppmD,SelParIdx);
		OutlierRem ot;
		if(ot.ORemMixModel(ppmD))
		{
			crt.cal_mean=ot.m_new;
			crt.cal_std=ot.std_new;
		}
		else
		{
			crt.cal_mean=ot.Mean(ppmD);
			crt.cal_std=ot.Std(crt.cal_mean,ppmD);
		}
		crt=PIso;
		bReturn.push_back(crt);	
	}
	else	
	{
		crt.charge=MS2Scans[MS2].ch;
		crt.cal_mean=MS2Scans[MS2].pmz-(ppmHV+ppmLV)/2;
		crt.iso_goodness=0.0;	
		crt.iso_num=0;
		crt.rel_int=0;
		crt.abs_int=0;
		crt.cal_std=(ppmHV-ppmLV)/6;
		bReturn.push_back(crt);
	}

	return true;
}


bool RawData::SCalibrateC_Lin(size_t MS2,vector<CalReturn>& bReturn)
{	
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;	
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;
	int MS1Scan=tmpIDX;	
	vector<CalData> PIso;
	//at first, we should find all the possible parent ions
	char buf[256];
	string stmp="Begin to Re-determinate the parent ions:\n";
	DInfBuf.push_back(stmp);

	clock_t start,finish; 
	start=clock();
	MSList[MS1Scan].ExtendFindISO(PIso,MS2Scans[MS2].pmz,MS2Scans[MS2].IW);
	finish=clock();
	sprintf(buf,"End to Re-determinate the parent ions:%lf\n",(double)(finish-start)/CLOCKS_PER_SEC);
	stmp=buf;
	DInfBuf.push_back(stmp);

	
#ifdef CASE_TEST
	int pionnum=PIso.size();
	fprintf(ParentIons,"%d\t%d\t%d\t%lf\n",pionnum,MS2Scans[MS2].scan,MS2Scans[MS2].ch,MS2Scans[MS2].pmz);	
#endif

	size_t pPNum=PIso.size();
	if(pPNum<=0) return false;
	vector<double>ppmD;
	CalReturn crt;
	for(size_t i=0;i<pPNum;i++)
	{
		//PIso[i].InstrumentType=instrument_T;
		LModel.Predict(PIso[i],ppmD,SelParIdx);
		OutlierRem ot;
		if(ot.ORemMixModel(ppmD))
		{
			crt.cal_mean=ot.m_new;
			crt.cal_std=ot.std_new;
		}
		else
		{
			crt.cal_mean=ot.Mean(ppmD);
			crt.cal_std=ot.Std(crt.cal_mean,ppmD);
		}
		crt=PIso[i];
		bReturn.push_back(crt);			
		ppmD.clear();
	}
	return true;
}

bool RawData::SCalibrateS_Lin(size_t MS2,vector<CalReturn>& bReturn)
{	
	bReturn.clear();
	if(MS2>=MS2Scans.size()) return false;	
	int tmpIDX=MS2Scans[MS2].MS1Scan;
	if(tmpIDX==-1) return false;
	int MS1Scan=tmpIDX;	
	CalData PIso;
	
	CalReturn crt;	
	if(MSList[MS1Scan].FindISO(PIso,MS2Scans[MS2].ch,MS2Scans[MS2].pmz))
	{
		vector<double>ppmD;		
		CalModel->Predict(PIso,ppmD,SelParIdx);
		OutlierRem ot;
		if(ot.ORemMixModel(ppmD))
		{
			crt.cal_mean=ot.m_new;
			crt.cal_std=ot.std_new;
		}
		else
		{
			crt.cal_mean=ot.Mean(ppmD);
			crt.cal_std=ot.Std(crt.cal_mean,ppmD);
		}
		crt=PIso;
		bReturn.push_back(crt);	
	}
	else	
	{
		crt.charge=MS2Scans[MS2].ch;
		crt.cal_mean=MS2Scans[MS2].pmz-(ppmHV+ppmLV)/2;
		crt.iso_goodness=0.0;	
		crt.iso_num=0;
		crt.rel_int=0;
		crt.abs_int=0;
		crt.cal_std=(ppmHV-ppmLV)/6;
		bReturn.push_back(crt);
	}
	return true;
}


double RawData::ECalibrate(size_t MS2)
{
	if(MS2>=MS2Scans.size()) return false;	
	int MS1Scan=MS2Scans[MS2].MS1Scan;
	if(MS1Scan==-1) return false;	
	CalData PIso;
	if(MSList[MS1Scan].FindISO(PIso,MS2Scans[MS2].ch,MS2Scans[MS2].pmz))
	{
		//PIso.InstrumentType=instrument_T;
		return CalModel->Predict(PIso,SelParIdx);
	}
	else return MS2Scans[MS2].pmz-(ppmHV+ppmLV)/2;
}

double RawData::ECalibrate_Lin(size_t MS2)
{
	if(MS2>=MS2Scans.size()) return false;		
	int MS1Scan=MS2Scans[MS2].MS1Scan;
	if(MS1Scan==-1) return false;	
	CalData PIso;
	if(MSList[MS1Scan].FindISO(PIso,MS2Scans[MS2].ch,MS2Scans[MS2].pmz))
	{
		//PIso.InstrumentType=instrument_T;
		return LModel.Predict(PIso,SelParIdx);		
	}
	else return MS2Scans[MS2].pmz-(ppmHV+ppmLV)/2;
}

void RawData::GeneratePurF(char tmpPureRawName[],char raw_names[])
{
	int idx=strlen(raw_names)-1;
	while(idx>=0)
	{
		if(raw_names[idx]=='\\') break;
		idx--;
	}
	idx++;
	if(idx<strlen(raw_names))strcpy(tmpPureRawName,raw_names+idx);
	else strcpy(tmpPureRawName,raw_names);

	//remove the .raw of .RAW
	idx=strlen(tmpPureRawName);
	char *pstr;
	if(idx>3)
	{
		pstr=tmpPureRawName+idx-4;
		if(strcmp(pstr,".raw")==0||strcmp(pstr,".RAW")==0)
		{
			pstr[0]='\0';
		}
	}	
}

bool RawData::IsParentExist(vector<CalReturn> &bReturn,double pmz)
{	
	size_t i,sg=bReturn.size();
	for(i=0;i<sg;i++)
	{
		if(bReturn[i].charge==CH_UNKNOWN) continue;
		double dmet=fabs(bReturn[i].cal_mean-pmz)*1e6/bReturn[i].cal_mean;	
		if(dmet<MET) return true;	
	}
	return false;
}


bool RawData::OutputMS2(FILE *MGFfp,FILE *IonsFile,size_t MS2,vector<CalReturn> &bReturn,char PureRawName[])
{	
	//string rawfile=raw_name;	
	//if(TIF==NULL)
	//{
	//	printf("the raw data can not be acessed.\n" );
	//	return false;
	//}
	//clock_t start,finish;
	string stmp;
	char buf[256];	
	sprintf(buf,"Begin to output MS2:%d.\n",MS2Scans[MS2].scan);
	stmp=buf;
	DInfBuf.push_back(stmp);
	//start=clock();
	Scan *ST;
	long scannum=MS2Scans[MS2].scan;
	ST=TIF->getScan(scannum);	
	size_t i,sg=bReturn.size();
	double MostPoss=0;
	double dm=0.5;

	bool IsPmzExist=IsParentExist(bReturn,MS2Scans[MS2].pmz);
	if(!IsPmzExist)//如果列表中没有有效的原母离子信息，直接输出一个
	{
		fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d.%d\n",PureRawName,ST->curScanNum,ST->curScanNum,MS2Scans[MS2].ch);	
		//bReturn[i].spectrum_num
		fprintf(MGFfp,"PEPMASS=%lf %lf\n",MS2Scans[MS2].pmz,ST->precursorIntensity_);//maybe wrong
		fprintf(MGFfp,"CHARGE=%d+\n",MS2Scans[MS2].ch);
		if(ms2format==MS2EMGF)
		{
			fprintf(MGFfp,"STOL=%lf\n",MET);
			fprintf(MGFfp,"MSNUM=0\n");
			//fprintf(MGFfp,"TOLU =ppm\n");
		}
		fprintf(MGFfp,"SCANS=%d\n",ST->curScanNum);
		fprintf(MGFfp,"RTINSECONDS=%lf\n",ST->retentionTimeInSec_);	
		ST->OutPutData(MGFfp);
		fprintf(MGFfp,"END IONS\n");
	}

	for(i=0;i<sg;i++)
	{
		if(bReturn[i].charge==CH_UNKNOWN) continue;//不输出未知电荷的情况；		
		//we will process this problem in the future
		//output
		double METInt;
		//if(IsMETModel)	METInt=3*METModel[1]*exp(bReturn[i].abs_int*METModel[3])+METModel[2];
		if(IsMETModel) METInt=3*(METModel[1]/sqrt(bReturn[i].spectrum_num*1.0)+METModel[2]);
		else  METInt=bReturn[i].cal_std*3e6/bReturn[i].cal_mean;
		if(METInt<1) ppbLevel++;
		//output ions specfic MET and other information
		fprintf(IonsFile,"%d\t%d\t%d\t%lf\t%lf\t%lf\n",bReturn[i].charge,ST->curScanNum,bReturn[i].spectrum_num,MS2Scans[MS2].pmz,bReturn[i].cal_mean,METInt);
		sprintf(buf,"Calibrate infor:MS_num=%d std=%lf scan=%d charge=%d\n",bReturn[i].spectrum_num,bReturn[i].cal_std,ST->curScanNum,bReturn[i].charge);
		stmp=buf;
		DInfBuf.push_back(stmp);
		sprintf(buf,"origianl pmz=%lf, new pmz=%lf\n",MS2Scans[MS2].pmz,bReturn[i].cal_mean);
		stmp=buf;
		DInfBuf.push_back(stmp);
		//fprintf(MGFfp,"PMASS=%lf\n",m*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
		if(MGFfp!=NULL)
		{
			//remedy for those unknown charges, we don not konw how the original
			//charge was assigned, may be according to the MS2 spectrum itself
		/*	if(bReturn[i].charge==CH_UNKNOWN)
			{
				double tdm=fabs(bReturn[i].cal_mean-ST->precursorMZ_)/ST->precursorMZ_;
				if(tdm<MET) bReturn[i].charge=ST->precursorCharge_;
			}*/
			fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d.%d\n",PureRawName,ST->curScanNum,ST->curScanNum,bReturn[i].charge);	
			//bReturn[i].spectrum_num
			fprintf(MGFfp,"PEPMASS=%lf %lf\n",bReturn[i].cal_mean,bReturn[i].abs_int);//maybe wrong
			fprintf(MGFfp,"CHARGE=%d+\n",bReturn[i].charge);
			if(ms2format==MS2EMGF)
			{
				fprintf(MGFfp,"STOL=%lf\n",METInt);
				fprintf(MGFfp,"MSNUM=%d\n",bReturn[i].spectrum_num);
				//fprintf(MGFfp,"TOLU =ppm\n");
			}
			fprintf(MGFfp,"SCANS=%d\n",ST->curScanNum);
			fprintf(MGFfp,"RTINSECONDS=%lf\n",ST->retentionTimeInSec_);	
			ST->OutPutData(MGFfp);
			fprintf(MGFfp,"END IONS\n");
		}

		//为了支持mzXML输出，得到和原来给出的mz最接近的一个输出
		double tdm=fabs(bReturn[i].cal_mean-ST->precursorMZ_);
		if(tdm<dm) 
		{
			dm=tdm;
			MostPoss=bReturn[i].cal_mean;
		}
	}
	delete ST;	

	if(MostPoss>10.0)
	{
		dm/=MostPoss;
		dm*=1e6;
		if(dm<MET)
		{
			TIF->scanmap.push_back(MS2Scans[MS2].scan);
			TIF->CalPmz.push_back(MostPoss);
		}
	}
	//finish=clock(); 	
	//double duration =(double)(finish-start)/CLOCKS_PER_SEC;
	//sprintf(buf,"Finished, the time consume: %lf.\n",duration);	
	//stmp=buf;
	//DInfBuf.push_back(stmp);
	return true;
}

void RawData::EnableOutput(char *fname)
{
	fp_cal_data=fopen(fname,"w");
	if(fp_cal_data!=NULL) IsOutput=true;
	else IsOutput=false;
}

void RawData::CloseOutput()
{
	if(IsOutput)
	{
		if(fp_cal_data!=NULL) fclose(fp_cal_data);
		IsOutput=false;
	}
}


int RawData::Loaddata(char *sraw)
{	
	strcpy(raw_name,sraw);//record the raw file name at first;
	string rawfile=sraw;
	if(TIF!=NULL) 
	{
		delete TIF;
		TIF=NULL;
	}
	TIF=new ThermoInterface;
	if (!TIF->initInterface())
	{
		printf("unable to interface with Thermo library" );
		delete TIF;
		TIF=NULL;
		return 0;
	}
	TIF->setCentroiding(true);
	if(!TIF->setInputFile(rawfile))
	{
		delete TIF;
		TIF=NULL;
		printf("raw file initial failure: %s.\n",sraw);	
		return 0;
	}
	//get the instrument type information
	instrument_T=TIF->instrumentInfo_.instrumentModel_;
	//suvery the header
	PreSuvey();
	//initial
	int CT[2];
	TIF->countMS(CT);
	if(CT[0]<=0||CT[1]<=0)
	{
		Error_code|=ERROR_RAWFAILUR;
		return 0;
	}
	MS2Scans.resize(CT[1]);
	MSList.resize(CT[0]);
	ScanMap.resize(CT[0]);
	size_t MS1Idx=0;
	size_t MS2Idx=0;
	//
	//
	clock_t start,finish; 
	printf("Finish initial the data file read interface,begin to process and count the time.\n");
	start=clock();
	Scan *ST;
	ST=TIF->getScan();	
	int bRetrun=0;	
	while(ST!=NULL)
	{
		if(TIF->MSLevel()==1)
		{
			MSList[MS1Idx].Convert(ST,ColParIdx,TailerIdx);
			MSList[MS1Idx].Pre_dm=MET;	
			ScanMap[MS1Idx]=ST->curScanNum;	
			MS1Idx++;				
		}
		else if(TIF->MSLevel()==2)
		{
			MS2Scans[MS2Idx].MS1Scan=MS1Idx-1;
			MS2Scans[MS2Idx].pmz=ST->precursorMZ_;
			MS2Scans[MS2Idx].IW=ST->IsolationWidth_;
			MS2Scans[MS2Idx].ch=ST->precursorCharge_;
			MS2Scans[MS2Idx].scan=ST->curScanNum;
			MS2Idx++;
			//debug
		/*	if(MS2Scan.scan==1914)
			{
				int kk=0;
				kk++;
			}*/
			//end						
		}	
		delete ST;
		ST=TIF->getScan();	
		bRetrun++;		
	}
	finish=clock(); 	
	double duration =(double)(finish-start)/CLOCKS_PER_SEC;
	printf("Finished, the time consume: %lf.\n",duration);		
	return bRetrun;
}

double RawData::GetSampleDm(double mass)
{
	if(instrument_T==LTQ_FT||instrument_T==LTQ_FT_ULTRA) return 6.06e-9*mass*mass-7.063e-11*mass+4.085e-8;//for FT
	else if(instrument_T==LTQ_ORBITRAP||instrument_T==LTQ_ORBITRAP_DISCOVERY||instrument_T==LTQ_ORBITRAP_XL||instrument_T==LTQ_ORBITRAP_VELOS)
		return 1.12e-009*mass*mass+2.204e-006*mass-0.0003366;//for Orbitrap
	else return 0;
}

void RawData::SortMerge(vector<mPeak> &tmpPKL)
{	
	//sort at first
	size_t i,k,sg=tmpPKL.size();
	mPeak mt;
	for(i=0;i<sg;i++)
	{
		for(k=i+1;k<sg;k++)
		{
			if(tmpPKL[i].dMass>tmpPKL[k].dMass) 
			{
				mt=tmpPKL[i];
				tmpPKL[i]=tmpPKL[k];
				tmpPKL[k]=mt;
			}
		}
	}
	//merge the neibour peaks
	k=0;
	vector<mPeak> Merg_PKL;
	while(k<sg)
	{	
		double dm_lim=GetSampleDm(tmpPKL[k].dMass)*5;//half peaks width
		for(i=k+1;i<sg;i++)
		{
			double dm=tmpPKL[i].dMass-tmpPKL[k].dMass;
			if(dm>dm_lim) break;
		}
		mt.dMass=0;
		mt.dIntensity=0;
		for(;k<i;k++)
		{
			mt.dMass+=tmpPKL[k].dMass*tmpPKL[k].dIntensity;
			mt.dIntensity+=tmpPKL[k].dIntensity;
		}
		if(mt.dIntensity>1e-2)mt.dMass/=mt.dIntensity;
		Merg_PKL.push_back(mt);		
	}
	tmpPKL.clear();
	sg=Merg_PKL.size();
	for(i=0;i<sg;i++) tmpPKL.push_back(Merg_PKL[i]);
}

void RawData::GetPossChFromPmz(double pmz,int scan,vector<int> &CHL)
{
	//determine the charge by combine MS within CH_RT_MIN, CH_RT_MAX
	vector<mPeak> Mix_PKL;
	int i,IDX=seekScan(scan,0,ScanMap.size());
	double min_mz=pmz-0.01;
	double max_mz=pmz+MAX_ISO*1.01;
	double Mid_RT=MSList[IDX].RT;
	for(i=IDX;i>=0;i--)
	{
		if(MSList[i].RT<CH_RT_MIN+Mid_RT) break;
		MSList[i].RetrivalPeaks(min_mz,max_mz,Mix_PKL);
	}
	int num_scan=MSList.size();

	for(i=IDX+1;i<num_scan;i++)
	{
		if(MSList[i].RT>CH_RT_MAX+Mid_RT) break;
		MSList[i].RetrivalPeaks(min_mz,max_mz,Mix_PKL);
	}
	SortMerge(Mix_PKL);
	//remove the irrelevant peaks
	size_t sg=Mix_PKL.size();
	CHL.clear();
	if(sg<=0) return;
	IsoGroup IG;
	sPeak pt=Mix_PKL[0];
	IG.Initial(pt);
	for(i=1;i<sg;i++) 
	{
		pt=Mix_PKL[i];
		IG.AddPeaks(pt,MET);
	}
	IG.GetCharge(CHL);//because the merge spectrum, 
	//the noise disterbut, can not perform overfit, just get the possible charges
}


bool RawData::Calibrate(char *rawfile,vector<Pre_calData> &IDCalData,char *mgfname)
{
	if(Error_code&ERROR_NOS_RFORMAT) return false;
	FILE *MGFfp,*IonsFile;
	if(ms2format!=MS2NON)
	{
		MGFfp=fopen(mgfname,"w");
		if(MGFfp==NULL)
		{
			//printf("can not create the output mgf file.\n");
			Error_code|=ERROR_COUPUTF;
			return false;
		}
	}
	else MGFfp=NULL;

	char tmpPureRawName[MAX_PATH];
	GeneratePurF(tmpPureRawName,raw_name);

	int slen=strlen(mgfname);
	mgfname[slen-1]='n';
	mgfname[slen-2]='o';
	mgfname[slen-3]='i';
	IonsFile=fopen(mgfname,"w");
	if(IonsFile==NULL)
	{
		fclose(MGFfp);		
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}

	if(Loaddata(rawfile)<=0)
	{
		//printf("can not load the raw data from file: %s.\n",rawfile);
		Error_code|=ERROR_RAWFAILUR;
		if(MGFfp!=NULL) fclose(MGFfp);
		return false;
	}	
	vector<CalData> XICList;	
	int i,pnum=IDCalData.size();	
	modeTotal=pnum;
	for(i=0;i<pnum;i++)
	{
		GetXICData(XICList,IDCalData[i].mze,IDCalData[i].charge,IDCalData[i].scan,IDCalData[i].mzt);
	}
	GetPreHis(XICList);	
	if(!ModelBuilding(XICList))
	{
		//printf("failure to build the calibration model.\n");
		Error_code|=ERROR_MDFAILUR;
		if(MGFfp!=NULL) fclose(MGFfp);
		return false;
	}
	GetCalHis(XICList);
	pnum=MS2Scans.size();
	vector<CalReturn> bReturn;
	CalTotal=0;
	for(i=0;i<pnum;i++)
	{
		if(XICCalibrate(i,bReturn))
		{
			OutputMS2(MGFfp,IonsFile,i,bReturn,tmpPureRawName);
			CalTotal+=bReturn.size();
		}		
	}
	if(MGFfp!=NULL) fclose(MGFfp);  
	fclose(IonsFile);
	return true;
}

//you must load data at first to use this function
bool RawData::Calibrate(char *mgfname)
{
	if(Error_code&ERROR_NOS_RFORMAT) return false;
	FILE *MGFfp,*IonsFile;
	if(ms2format!=MS2NON)
	{
		MGFfp=fopen(mgfname,"w");
		if(MGFfp==NULL)
		{
			//printf("can not create the output mgf file.\n");
			Error_code|=ERROR_COUPUTF;
			return false;
		}	
	}
	else MGFfp=NULL;

	char tmpPureRawName[MAX_PATH];
	GeneratePurF(tmpPureRawName,raw_name);

	int slen=strlen(mgfname);
	mgfname[slen-1]='n';
	mgfname[slen-2]='o';
	mgfname[slen-3]='i';
	IonsFile=fopen(mgfname,"w");
	if(IonsFile==NULL)
	{
		if(MGFfp!=NULL) fclose(MGFfp);		
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}

	size_t i, pnum=MS2Scans.size();
	vector<CalReturn> bReturn;
	CalTotal=0;
	for(i=0;i<pnum;i++)
	{
		//to debug the infer-cycle problem
		//if(MS2Scans[i].scan==298)
		//{
		//	int kk=0;
		//	kk++;
		//}

		if(XICCalibrate(i,bReturn))
		{
			OutputMS2(MGFfp,IonsFile,i,bReturn,tmpPureRawName);
			 CalTotal++;
		}		
	}
	if(MGFfp!=NULL)  fclose(MGFfp); 
	fclose(IonsFile);
	return true;
}
//you must load data at first to use this function
bool RawData::Calibrate_SVM(char *mgfname)
{
#ifdef CASE_TEST
	//FILE *ParentIons;
	char pfname[MAX_PATH];
	strcpy(pfname,mgfname);
	int slen=strlen(mgfname);
	pfname[slen-3]='t';
	pfname[slen-2]='o';
	pfname[slen-1]='n';
	ParentIons=fopen(pfname,"w");
	if(ParentIons==NULL)
	{
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}	
#endif

	//GPOS(&ConselRow,&ConselCol);
	DInfBuf.clear();
	bool bReturn;
	//if(IS_DYA_XIC) DFT.Initial(DYA_XIC_SW,DYA_XIC_SM,2);
	if(XICCALIBRATE) bReturn= Calibrate_Complex(mgfname);
	else bReturn= Calibrate_Simple(mgfname);

#ifdef CASE_TEST
	//FILE *ParentIons;
	fclose(ParentIons);	
#endif
	return bReturn;
}

//you must load data at first to use this function
bool RawData::Calibrate_Non(char *mgfname)
{
	if(Error_code&ERROR_NOS_RFORMAT) return false;

	//if(IS_DYA_XIC) DFT.Initial(DYA_XIC_SW,DYA_XIC_SM,2);

	FILE *MGFfp,*IonsFile;
	if(ms2format!=MS2NON)
	{
		MGFfp=fopen(mgfname,"w");
		if(MGFfp==NULL)
		{
			//printf("can not create the output mgf file.\n");
			Error_code|=ERROR_COUPUTF;
			return false;
		}
	}
	else MGFfp=NULL;

	char tmpPureRawName[MAX_PATH];
	GeneratePurF(tmpPureRawName,raw_name);

	int slen=strlen(mgfname);
	mgfname[slen-1]='n';
	mgfname[slen-2]='o';
	mgfname[slen-3]='i';
	IonsFile=fopen(mgfname,"w");
	if(IonsFile==NULL)
	{
		if(MGFfp!=NULL) fclose(MGFfp);	
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}

	size_t i,pnum;
	pnum=MSList.size();
	double tmpdm;
	if(PRE_DM<1e-2) tmpdm=(ppmHV-ppmLV)/SQRTTWO;
	else tmpdm=PRE_DM;
	for(i=0;i<pnum;i++) MSList[i].Pre_dm=tmpdm;

	pnum=MS2Scans.size();
	vector<CalReturn> bReturn;
	CalTotal=0;
	for(i=0;i<pnum;i++)
	{		
		if(XICCalibrate_Non(i,bReturn))
		{
			OutputMS2(MGFfp,IonsFile,i,bReturn,tmpPureRawName);
			CalTotal+=bReturn.size();
		}
		OutPutInf();

	}
	if(MGFfp!=NULL) fclose(MGFfp);  
	fclose(IonsFile);
	return true;
}

//you must load data at first to use this function
bool RawData::Calibrate_Lin(char *mgfname)
{
#ifdef CASE_TEST
	//FILE *ParentIons;
	char pfname[MAX_PATH];
	strcpy(pfname,mgfname);
	int slen=strlen(mgfname);
	pfname[slen-3]='t';
	pfname[slen-2]='o';
	pfname[slen-1]='n';
	ParentIons=fopen(pfname,"w");
	if(ParentIons==NULL)
	{
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}	
#endif

	//GPOS(&ConselRow,&ConselCol);
	DInfBuf.clear();
	bool bReturn;
	//if(IS_DYA_XIC) DFT.Initial(DYA_XIC_SW,DYA_XIC_SM,2);
	if(XICCALIBRATE) bReturn= Calibrate_Com_Lin(mgfname);
	else bReturn= Calibrate_Sim_Lin(mgfname);

#ifdef CASE_TEST
	//FILE *ParentIons;
	fclose(ParentIons);	
#endif
	return bReturn;
}

//you must load data at first to use this function
bool RawData::Calibrate_Complex(char *mgfname)
{
	if(Error_code&ERROR_NOS_RFORMAT) return false;
	FILE *MGFfp,*IonsFile;
	if(ms2format!=MS2NON)
	{
		MGFfp=fopen(mgfname,"w");
		if(MGFfp==NULL)
		{
			//printf("can not create the output mgf file.\n");
			Error_code|=ERROR_COUPUTF;
			return false;
		}	
	}
	else MGFfp=NULL;
	
	char tmpPureRawName[MAX_PATH];
	GeneratePurF(tmpPureRawName,raw_name);

	int slen=strlen(mgfname);
	mgfname[slen-1]='n';
	mgfname[slen-2]='o';
	mgfname[slen-3]='i';
	IonsFile=fopen(mgfname,"w");
	if(IonsFile==NULL)
	{
		if(MGFfp!=NULL)  fclose(MGFfp);		
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}

	size_t i,pnum;
	pnum=MSList.size();
	double tmpdm;
	if(PRE_DM<1e-2) tmpdm=(ppmHV-ppmLV)/SQRTTWO;
	else tmpdm=PRE_DM;
	for(i=0;i<pnum;i++) MSList[i].Pre_dm=tmpdm;

	pnum=MS2Scans.size();
	vector<CalReturn> bReturn;
	CalTotal=0;
	for(i=0;i<pnum;i++)
	{
		//to debug
		//if(MS2Scans[i].scan==3391)//1914
		//{
		//	int kk=0;
		//	kk++;
		//}
		////		
		if(XICCalibrate(i,bReturn))
		{
			OutputMS2(MGFfp,IonsFile,i,bReturn,tmpPureRawName);
			CalTotal+=bReturn.size();

		}
		OutPutInf();
	}
	if(MGFfp!=NULL)  fclose(MGFfp);  
	fclose(IonsFile);
	return true;
}

//you must load data at first to use this function
bool RawData::Calibrate_Simple(char *mgfname)
{
	if(Error_code&ERROR_NOS_RFORMAT) return false;
	FILE *MGFfp,*IonsFile;
	if(ms2format!=MS2NON)
	{
		MGFfp=fopen(mgfname,"w");
		
		if(MGFfp==NULL)
		{
			//printf("can not create the output mgf file.\n");
			Error_code|=ERROR_COUPUTF;
			return false;
		}	
	}
	else MGFfp=NULL;

	char tmpPureRawName[MAX_PATH];
	GeneratePurF(tmpPureRawName,raw_name);

	int slen=strlen(mgfname);
	mgfname[slen-1]='n';
	mgfname[slen-2]='o';
	mgfname[slen-3]='i';
	IonsFile=fopen(mgfname,"w");
	if(IonsFile==NULL)
	{
		if(MGFfp!=NULL)  fclose(MGFfp);		
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}

	size_t i,pnum;
	pnum=MSList.size();
	double tmpdm;
	if(PRE_DM<1e-2) tmpdm=(ppmHV-ppmLV)/SQRTTWO;
	else tmpdm=PRE_DM;
	for(i=0;i<pnum;i++) MSList[i].Pre_dm=tmpdm;

	pnum=MS2Scans.size();
	vector<CalReturn> bReturn;
	CalTotal=0;
	for(i=0;i<pnum;i++)
	{
		//to debug
		//if(MS2Scans[i].scan==3391)//1914
		//{
		//	int kk=0;
		//	kk++;
		//}
		////		
		if(SCalibrate(i,bReturn))
		{
			OutputMS2(MGFfp,IonsFile,i,bReturn,tmpPureRawName);
			CalTotal+=bReturn.size();
		}
		OutPutInf();

	}
	if(MGFfp!=NULL) fclose(MGFfp);  
	fclose(IonsFile);
	return true;
}

//you must load data at first to use this function
bool RawData::Calibrate_Com_Lin(char *mgfname)
{
	if(Error_code&ERROR_NOS_RFORMAT) return false;
	FILE *MGFfp,*IonsFile;
	if(ms2format!=MS2NON)
	{
		MGFfp=fopen(mgfname,"w");
		if(MGFfp==NULL)
		{
			//printf("can not create the output mgf file.\n");
			Error_code|=ERROR_COUPUTF;
			return false;
		}	
	}
	else MGFfp=NULL;

	char tmpPureRawName[MAX_PATH];
	GeneratePurF(tmpPureRawName,raw_name);

	int slen=strlen(mgfname);
	mgfname[slen-1]='n';
	mgfname[slen-2]='o';
	mgfname[slen-3]='i';
	IonsFile=fopen(mgfname,"w");
	if(IonsFile==NULL)
	{
		if(MGFfp!=NULL) fclose(MGFfp);		
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}


	size_t i,pnum;
	pnum=MSList.size();
	double tmpdm;
	if(PRE_DM<1e-2) tmpdm=(ppmHV-ppmLV)/SQRTTWO;
	else tmpdm=PRE_DM;
	for(i=0;i<pnum;i++) MSList[i].Pre_dm=tmpdm;

	pnum=MS2Scans.size();
	vector<CalReturn> bReturn;
	CalTotal=0;
	for(i=0;i<pnum;i++)
	{
		//to debug
		//if(MS2Scans[i].scan==3391)//1914
		//{
		//	int kk=0;
		//	kk++;
		//}
		////		
		if(XICCalibrate_Lin(i,bReturn))
		{
			OutputMS2(MGFfp,IonsFile,i,bReturn,tmpPureRawName);
			CalTotal+=bReturn.size();

		}
		OutPutInf();
	}
	if(MGFfp!=NULL) fclose(MGFfp);  
	fclose(IonsFile);
	return true;
}

//you must load data at first to use this function
bool RawData::Calibrate_Sim_Lin(char *mgfname)
{
	if(Error_code&ERROR_NOS_RFORMAT) return false;
	FILE *MGFfp,*IonsFile;
	if(ms2format!=MS2NON)
	{
		MGFfp=fopen(mgfname,"w");
		if(MGFfp==NULL)
		{
			//printf("can not create the output mgf file.\n");
			Error_code|=ERROR_COUPUTF;
			return false;
		}	
	}
	else MGFfp=NULL;

	char tmpPureRawName[MAX_PATH];
	GeneratePurF(tmpPureRawName,raw_name);


	int slen=strlen(mgfname);
	mgfname[slen-1]='n';
	mgfname[slen-2]='o';
	mgfname[slen-3]='i';
	IonsFile=fopen(mgfname,"w");
	if(IonsFile==NULL)
	{
		if(MGFfp!=NULL) fclose(MGFfp);		
		//printf("can not create the output mgf file.\n");
		Error_code|=ERROR_COUPUTF;
		return false;
	}

	size_t i,pnum;
	pnum=MSList.size();
	double tmpdm;
	if(PRE_DM<1e-2) tmpdm=(ppmHV-ppmLV)/SQRTTWO;
	else tmpdm=PRE_DM;
	for(i=0;i<pnum;i++) MSList[i].Pre_dm=tmpdm;

	pnum=MS2Scans.size();
	vector<CalReturn> bReturn;
	CalTotal=0;
	for(i=0;i<pnum;i++)
	{
		//to debug
		//if(MS2Scans[i].scan==3391)//1914
		//{
		//	int kk=0;
		//	kk++;
		//}
		////		
		if(SCalibrate_Lin(i,bReturn))
		{
			OutputMS2(MGFfp,IonsFile,i,bReturn,tmpPureRawName);
			CalTotal+=bReturn.size();
		}
		OutPutInf();
	}
	if(MGFfp!=NULL) fclose(MGFfp); 
	fclose(IonsFile);
	return true;
}

int RawData::GetErrorCode()
{
	return Error_code;
}

void RawData::SetErrorCode(int errors)
{
	Error_code|=errors;
}
//需要更新中，无模型校正也需要建模，以估计误差范围
void RawData::GenReport(CReport *rep)
{
	if(Error_code&ERROR_RAWFAILUR||Error_code&ERROR_NOTLOAD) return;//根本没有进行建模，不需要产生报告了。
	rep->count=CalTotal;
	rep->mcount=modeTotal;
	rep->mtA=mtA;
	rep->mtB=mtB;
	rep->Error_code=Error_code;
	rep->IsMETModel=IsMETModel;
	if(ModelType==MODEL_SVM)
	{
		rep->Model_mean=CalModel->Res_mean;
		rep->Model_std=CalModel->Res_std;
		rep->Model_gd=CalModel->gd;
	}
	else if(ModelType==MODEL_LIN)
	{
		rep->Model_mean=LModel.stats[0];
		rep->Model_std=LModel.stats[1];
		rep->Model_gd=LModel.stats[2];
		for(int i=0;i<LModel.PARNUM;i++)
		{
			rep->LCoff.push_back(LModel.COFF[i]);
		}
	}
	else
	{
		rep->Model_mean=0;
		rep->Model_std=0;
		rep->Model_gd=0;
	}

	size_t i,sg=SelParIdx.size();
	rep->SelParName.clear();
	rep->SelParName.push_back("mz experiment");
	//if(ModelType==MODEL_SVM) rep->SelParName.push_back("mass experiment");
	rep->SelParName.push_back("parent ion intensity");
	rep->SelParName.push_back("total ion current");
	//if(ModelType==MODEL_SVM) rep->SelParName.push_back("isotopic peak number");	
	rep->SelParName.push_back("retention time");
	rep->SelParName.push_back("signal/noise");
		rep->SelParName.push_back("relative intensity");
	for(i=0;i<sg;i++)
	{
		rep->SelParName.push_back(ParNames[SelParIdx[i]]);
	}

	for(i=0;i<4;i++) rep->METModel[i]=METModel[i];

	rep->ppbLevel=ppbLevel;
}

long RawData::GetScanNum(CString sFName)
{
	int slen=sFName.GetLength()-1;
	CString item;
	while(slen>=0&&sFName[slen]!='.') slen--;
	slen--;
	while(slen>=0&&sFName[slen]!='.') slen--;
	slen--;
	int num=0;
	while(slen>=0&&sFName[slen]!='.') 
	{
		slen--;
		num++;
	}
	if(slen>=0) 
	{
		slen++;
		CString item=sFName.Mid(slen,num);
		return atoi(item);
	}
	return -1;
}

int RawData::CollectDataS(CString RPath,vector<Pre_calData> &IDCalData)
{	
	WIN32_FIND_DATA dat;
    HANDLE hand;
	CString tfname,sTemp;
	tfname=RPath+"\\*.out";
	hand=::FindFirstFile(tfname,&dat);//find the first *.out file
	if(hand==INVALID_HANDLE_VALUE) 
	{
		Error_code|=ERROR_RNOPEN;	
		return 0;
	}
	int count=0;
	Pre_calData pt;
	SRecord rec;
	OutRecord OutFile;
	tfname=RPath+'\\'+dat.cFileName;
	if(OutFile.ReadFromFile(tfname))
	{
		OutFile.GetAt(0,&rec);	
		rec.detCn=OutFile.GetdetCn();
		if(cFilter1.IsPassFilter(&rec,OutFile.Charge))
		{
			pt.mze=(OutFile.EMH+(OutFile.Charge-1)*1.007825f)/OutFile.Charge;
			pt.mzt=(rec.MH+(OutFile.Charge-1)*1.007825f)/OutFile.Charge;
			pt.charge=OutFile.Charge;	
			pt.scan=GetScanNum(tfname);	
			IDCalData.push_back(pt);
			count++;
		}		
	}

	while(::FindNextFile(hand,&dat))//explore all the file in this directory
	{
		tfname=RPath+'\\'+dat.cFileName;
		if(OutFile.ReadFromFile(tfname))
		{
			OutFile.GetAt(0,&rec);	
			rec.detCn=OutFile.GetdetCn();
			if(cFilter1.IsPassFilter(&rec,OutFile.Charge))
			{
				pt.mze=(OutFile.EMH+(OutFile.Charge-1)*1.007825f)/OutFile.Charge;
				pt.mzt=(rec.MH+(OutFile.Charge-1)*1.007825f)/OutFile.Charge;
				pt.charge=OutFile.Charge;	
				pt.scan=GetScanNum(tfname);	
				IDCalData.push_back(pt);
				count++;
			}	
		}
	}
	::FindClose(hand);//关闭查找文件需要的句柄			
	if(count==0) Error_code|=ERROR_NOTLOAD;	
	return count;
}

int RawData::CollectDataP(CString RFile,vector<Pre_calData> &IDCalData)
{
	PepXML PRT;
	char buf[1024];
	strcpy(buf,RFile.GetBuffer());
	RFile.ReleaseBuffer();	
	int count=PRT.LoadData(buf,IDCalData);
	if(count==0) Error_code|=ERROR_RNOPEN;	
	return count;	
}

int RawData::CollectDataX(CString RFile,vector<Pre_calData> &IDCalData)
{
	MassMatrixR MRF;
	char buf[1024];
	strcpy(buf,RFile.GetBuffer());
	RFile.ReleaseBuffer();	
	int count=MRF.LoadData(buf,IDCalData);
	if(count==0) Error_code|=ERROR_RNOPEN;	
	return count;	
}

int RawData::CollectDataM(CString RFile,vector<Pre_calData> &IDCalData)
{
	MResult MascotR;
	MascotR.cFilter=cFilter;
	if(!MascotR.LoadData(RFile)) 
	{ 
		Error_code|=ERROR_RNOPEN;	
		return 0;
	}
	int count=MascotR.GetCount();		
	iPeptide pt;
	Pre_calData dt;
	int num=0;
	for(int i=0;i<count;i++)
	{
		MascotR.GetAt(i,&pt);
		if(pt.ScanNumber<=0) continue;
		dt.mze=pt.Observed;	
		dt.charge=(int)((pt.ObsrvMass+MASS_H)/pt.Observed+0.3);
		dt.mzt=pt.CalcMass/(dt.charge)+MASS_H_ME;	
		dt.scan=pt.ScanNumber;	
		IDCalData.push_back(dt);
		num++;
	}
	if(num==0)	Error_code|=ERROR_NOTLOAD;
	return num;
}

void RawData::InitialFilter(Filter *dFilter,SFilter *dFilter1)
{
	cFilter=*dFilter;
	cFilter1=*dFilter1;
}

bool RawData::Initial(char *rawfname)
{
	if(Loaddata(rawfname)<=0)
	{
		//printf("can not load the raw data from file: %s.\n",rawfile);
		Error_code|=ERROR_RAWFAILUR;	
		return false;
	}
	mtA.rawname=rawfname;
	mtB.rawname=rawfname;
	ppbLevel=0;
	return true;

}

void RawData::OutputModel(char *svmmodel)
{
	CalModel->OutputModel(svmmodel);
}

void RawData::SetOutputRData(char *fname)
{
	strcpy(IntRModelData,fname);
	OIntR_Data=true;
}

//to clear for a new states
void RawData::clear()
{	
	if(IsOutput)
	{
		if(fp_cal_data!=NULL) fclose(fp_cal_data);
		fp_cal_data=NULL;
	}
	if(TIF!=NULL) 
	{
		delete TIF;	
		TIF=NULL;
	}

	if(CalModel!=NULL)
	{
		delete CalModel;
		CalModel=NULL;
	}

	raw_name[0]='\0';
	for(int i=0;i<MAX_ISO;i++) IsoDisT[i]=0;
	MET=10;
		
	instrument_T=INSTRUMENTMODEL_UNDEF;
	Error_code=0;
	CalTotal=0;
	modeTotal=0;
	ppmLV=-10;
	ppmHV=10;
	METModel[0]=0;
	METModel[1]=0;
	METModel[2]=0;
	METModel[3]=0;
	IsMETModel=false;

	IntRModelData[0]='\0';
	OIntR_Data=false;
}

SVMModel *RawData::GetSVMPointer()
{
	return CalModel;
}

void RawData::PreSuvey()
{
	Scan *ST;
	ST=TIF->getScan();
	while(ST!=NULL)
	{
		if(TIF->MSLevel()==1) break;
		delete ST;		
		ST=TIF->getScan();			
	}
	ColParIdx.clear();
	ParNames.clear();	
	//ParNames.push_back("Ion Injection Time");
	//ParNames.push_back("Elapsed Scan Time");
	size_t i,sg=ST->status_par_name.size();
	double tvalue;
	for(i=0;i<sg;i++)
	{
		if(sscanf(ST->status_par_value[i].c_str(),"%lf",&tvalue)==1)
		{
			ParNames.push_back(ST->status_par_name[i]);
			ColParIdx.push_back(i);
		}
	}
	sg=ST->tailer_par_name.size();
	for(i=0;i<sg;i++)
	{
		size_t IDX=ST->tailer_par_name[i].find("Ion Injection Time");
		size_t IDX1=ST->tailer_par_name[i].find("Elapsed Scan Time");
		if(IDX!=string::npos)TailerIdx[0]=IDX;
		if(IDX1!=string::npos)TailerIdx[1]=IDX1;
	}
	delete ST;
	TIF->SeekToScan(TIF->getFirstScanNumber());
	//cal and copy the mz interval model
	TIF->PreSuryMZINT();
	for(int i=0;i<4;i++) MzIntModel[i]=TIF->MzIntModel[i];	
}

void RawData::featureSel(vector<CalData> &XICList)
{
	if(!IsUseExtendPar)
	{
		SelParIdx.clear();
		return;
	}
	size_t TC=XICList.size();
	if(TC<=0) return;
	size_t Col=XICList[0].FTStatus.size();
	if(Col<=0) return;
	double *feav,*target,*Corr;
	feav=new double[TC];
	target=new double[TC];
	Corr=new double[Col];	
	size_t i,k;
	for(i=0;i<TC;i++) target[i]=(XICList[i].Isomz[0]-XICList[i].mzt)*1e6/XICList[i].mzt;
	for(i=0;i<Col;i++)
	{
		for(k=0;k<TC;k++) feav[k]=XICList[k].FTStatus[i];
		double m_f=gsl_stats_mean(feav,1,TC);
		double std_f=gsl_stats_sd_m(feav,1,TC,m_f);
		if(std_f<0.01*fabs(m_f)||std_f<1e-6) Corr[i]=0;
		else Corr[i]=fabs(gsl_stats_correlation(feav,1,target,1,TC));
	}
	double mc=gsl_stats_mean(Corr,1,Col);
	double sd=gsl_stats_sd_m(Corr,1,Col,mc);
	double CT=mc+sd;
	SelParIdx.clear();
	for(i=0;i<Col;i++) 
	{
		if(Corr[i]>CT) SelParIdx.push_back(i);
	}
	delete []Corr;
	delete []feav;
	delete []target;
}

void RawData::OutPutInf()
{
	//SPOS(ConselRow,0);
	size_t s=DInfBuf.size();
	for(size_t i=0;i<s;i++)
	{
		printf((TCHAR*)DInfBuf[i].c_str());
	}
	DInfBuf.clear();
}

double RawData::GetSmDm(double mass)
{
	if(MzIntModel[3]>0.95)
	{
		return MzIntModel[0]+MzIntModel[1]*mass+MzIntModel[2]*mass*mass;
	}
	else return GetSampleDm(mass);
}

bool RawData::SaveRawData(char *fname)
{
	clock_t   startclock,   endclock;  
	startclock=clock(); 
	FILE *fp;
	fp=fopen(fname,"wb");
	if(fp==NULL) return false;
	//写入raw文件名
	size_t i,size=strlen(raw_name);
	fwrite(&size,sizeof(size_t),1,fp);
	fwrite(raw_name,sizeof(char),size,fp);
	//写入仪器类型
	fwrite(&instrument_T,sizeof(MSInstrumentModelType),1,fp);
	//写入状态参数名称
	size=ParNames.size();
	fwrite(&size,sizeof(size_t),1,fp);
	size_t count;
	for(i=0;i<size;i++)
	{
		count=ParNames[i].size();
		fwrite(&count,sizeof(size_t),1,fp);
		fwrite(ParNames[i].c_str(),sizeof(char),count,fp);
	}
	//写入scanmap
	size=ScanMap.size();
	fwrite(&size,sizeof(size_t),1,fp);
	for(i=0;i<size;i++)
	{
		fwrite(&(ScanMap[i]),sizeof(long),1,fp);
	}
	//写入MS2header数据
	size=MS2Scans.size();
	fwrite(&size,sizeof(size_t),1,fp);
	for(i=0;i<size;i++)
	{
		MS2Scans[i].Write2File(fp);
	}
	//写入图谱数据
	size=MSList.size();	
	fwrite(&size,sizeof(size_t),1,fp);
	for(i=0;i<size;i++)
	{
		MSList[i].Write2File(fp);
	}
	fclose(fp);
	endclock=clock(); 
	printf("The time consuume is :%f\n",(float)(endclock-startclock)/CLOCKS_PER_SEC); 
	return true;
}

bool RawData::LoadRawData(char *fname)
{
	clock_t   startclock,   endclock;  
	startclock=clock();  	
	FILE *fp;
	fp=fopen(fname,"rb");
	if(fp==NULL) return false;
	size_t i,count=0;
	//读取文件名
	fread(&count,sizeof(size_t),1,fp);
	fread(raw_name,sizeof(char),count,fp);
	//读出仪器类型
	fread(&instrument_T,sizeof(MSInstrumentModelType),1,fp);
	//读出状态参数名
	fread(&count,sizeof(size_t),1,fp);	
	char buf[1024];
	string st;
	ParNames.clear();
	for(i=0;i<count;i++)
	{	
		size_t sg;
		fread(&sg,sizeof(size_t),1,fp);
		fread(buf,sizeof(char),sg,fp);
		buf[sg]='\0';
		st=buf;
		ParNames.push_back(st);
	}
	//读出scan map数据
	fread(&count,sizeof(size_t),1,fp);	
	ScanMap.clear();
	long tscan;
	for(i=0;i<count;i++)
	{
		fread(&tscan,sizeof(long),1,fp);
		ScanMap.push_back(tscan);
	}
	//读出MS2header数据
	fread(&count,sizeof(size_t),1,fp);
	MS2Header mt;
	MS2Scans.clear();
	for(i=0;i<count;i++)
	{
		mt.ReadFromFile(fp);
		MS2Scans.push_back(mt);
	}
	//读出图谱数据
	fread(&count,sizeof(size_t),1,fp);
	MSList.clear();	
	CalScan tMS1;
	for(i=0;i<count;i++)
	{
		tMS1.ReadFromFile(fp);
		MSList.push_back(tMS1);
	}
	endclock=clock(); 
	printf("The time consuume is :%f\n",(float)(endclock-startclock)/CLOCKS_PER_SEC); 
	return true;
}


void RawData::GetPurRawName(string &purname)
{	
	purname=raw_name;
	size_t found=purname.rfind('.');
	if(found!=string::npos)
	{
		purname=purname.substr(0,found);
	}	

	found=purname.rfind('\\');
	if(found!=string::npos)
	{
		purname=purname.substr(found+1);
	}
}

bool RawData::OutPutMS1(string outpath)
{	
	if(TIF==NULL)return false;	
	TIF->MType=ModelType;	
	TIF->IsCalibrate=true;
	TIF->ReInitial();
	if(ModelType==MODEL_SVM) TIF->CalModel=GetSVMPointer();
	else if(ModelType==MODEL_LIN)TIF->LM=LModel;
	else TIF->IsCalibrate=false;

	TIF->totalpmzscan=TIF->scanmap.size();
	TIF->totalpmzscan--;
	//TIF->IsUseCpmz=true;
	TIF->Convert(ColParIdx,SelParIdx,TailerIdx);
	ConverterArgs args;
	args.centroidScans=true;
	args.compressScans=false;
	args.gzipOutputFile = false;
	args.verbose=true;	
	args.forcePrecursorFromFilterLine=false;
	args.shotgunFragmentation = false;
	args.lockspray = false;
	string tmps;
	GetPurRawName(tmps);
	args.outputFileName=outpath;	
	args.outputFileName+="\\";
	args.outputFileName+=tmps;	
	if(ms1format==MS1_mzML)
	{
		args.mzXMLMode = false;
		args.mzMLMode = true;
		args.outputFileName += ".mzML";
	}
	else
	{
		args.mzXMLMode = true;
		args.mzMLMode = false;
		args.outputFileName += ".mzXML";
	}
	args.threshold = false;
	args.thresholdDiscard = false;	
	args.inclusiveCutoff=0;	
	args.inputFileName=raw_name;
	return  msconvert(args,*TIF);
}

void RawData::InitialDFTModel()
{
	DFT.Initial(DYA_XIC_SW,DYA_XIC_SM,2);
}