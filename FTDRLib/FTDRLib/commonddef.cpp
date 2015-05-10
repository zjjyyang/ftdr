#include "stdafx.h"
#include "commonddef.h"
#include "math.h"

int ConselRow=0;//保存当前调试窗口光标的行位置
int ConselCol=0;//保存当前调试窗口光标的列位置

double MINRT=-5.0;
double MAXRT=5.0;
double GD_CUT=0.2;//同位素峰的拟合优度门限
int MAX_CH=4;//最大的考虑电荷
//int IT_LIM=5;//XIC构建中，信号缺失的最大次数
int INT_TIME_MAX=3;
double RINT_CUT=0.001;//相对信号强度过滤门限
int ISO_CUT=1;//means no cut here同位素峰数目的限制
double MAX_INT_RT=0.2;
double MIN_SN=2;//信噪比门限，默认为2

double CH_RT_MIN=-5;//电荷确定中，搜索图谱叠加的区间向左扩展范围
double CH_RT_MAX=5;//电荷确定中，搜索图谱叠加的区间向右扩展范围

long MAX_SVM_TRAIN=6000;//SVM模型训练样本数量最大值

int MIN_H_ISO_NUM=2;//扩展查找母离子时，同位素峰最小值
double R_INT_CUT_H= 0.01;//扩展查找时，相对信号强度最小值
double GD_CUT_H= 0.2;//original is 0.2

int LLEXTEND=10;// 10
int RLEXTEND=10;// 10

double PRE_DM=0;//10ppm
double PRE_DM_10=15;
bool XICCALIBRATE=true;//true
bool PARIONRED=true;//true;
double SVMC=64;

int DYA_XIC_SW=21;
int DYA_XIC_SM=7;
bool IS_DYA_XIC=false;

ThreadParm::ThreadParm()
{
	IsEnd=false;
	ThreadID=0;
	hEventP=NULL;
	hEventE=NULL;
	RepList=NULL;
	hEventB=NULL;
	count=0;
	IsDeal=false;
	IsOpCData=false;	

	ms1format=MS1_NON;
	ms2format=MS2MGF;

	IsUseStatsPar=true;

	ModelType=MODEL_SVM;
}

ThreadParm::~ThreadParm()
{
	if(RepList!=NULL&&count>0) delete []RepList;
}

ThreadParm::ThreadParm(const ThreadParm &r)
{
	RFile=r.RFile;
	RawFile=r.RawFile;
	sOutPath=r.sOutPath;
	IsEnd=r.IsEnd;
	ThreadID=r.ThreadID;
	dFilter=r.dFilter;
	dFilter1=r.dFilter1;
	type=r.type;
	RepList=r.RepList;
	count=r.count;
	IsDeal=r.IsDeal;	
	IsOpCData=r.IsOpCData;
	IsOPSVMModel=r.IsOPSVMModel;
	//SVMModelFile=r.SVMModelFile;	
	ms1format=r.ms1format;
	ms2format=r.ms2format;

	IsUseStatsPar=r.IsUseStatsPar;
	ModelType=r.ModelType;
}
ThreadParm &ThreadParm::operator=(const ThreadParm &r)
{
	RFile=r.RFile;
	RawFile=r.RawFile;
	sOutPath=r.sOutPath;
	IsEnd=r.IsEnd;
	ThreadID=r.ThreadID;
	dFilter=r.dFilter;
	dFilter1=r.dFilter1;
	type=r.type;
	RepList=r.RepList;
	count=r.count;
	IsDeal=r.IsDeal;	
	IsOpCData=r.IsOpCData;
	IsOPSVMModel=r.IsOPSVMModel;
	//SVMModelFile=r.SVMModelFile;	
	ms1format=r.ms1format;
	ms2format=r.ms2format;

	IsUseStatsPar=r.IsUseStatsPar;
	ModelType=r.ModelType;
	return *this;
}
void ThreadParm::AllocRepBuff(int n)
{
	if(RepList!=NULL&&count>0) 
	{
		delete []RepList;
		count=0;
		RepList=NULL;
	}
	if(n>0)
	{
		RepList=new CReport[n];	
		count=n;
	}
}

void ThreadParm::FreeRepBuff()
{
	if(RepList!=NULL&&count>0) 
	{
		delete []RepList;
		count=0;
		RepList=NULL;
	}
}

bool ThreadParm::SetRepBuff(int idx,CReport *rt)
{
	if(count>idx&&idx>=0)
	{
		RepList[idx]=*rt;
		return true;
	}
	return false;
}

bool ThreadParm::GetRepBuff(int idx,CReport *rt)
{
	if(count>idx&&idx>=0)
	{
		*rt=RepList[idx];
		return true;
	}
	return false;
}

void ThreadParm::ResetRepState(int idx)
{
	if(count>idx&&idx>=0)
	{
		RepList[idx].IsRepValidate=false;
	}
}
///////////////////////////////////////////////////////////////////////
CalData::CalData()
{
	//all the initial values is 0
	goodness=0;
	IsoNum=0;
	mzt=0;
	scannum=0;
	tic=0;
	ch=0;
	RT=0;
	BaseInt=0;
	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		IsotopicE[i]=0;
		IsotopicT[i]=0;
		Isomz[i]=0;
	}
}

CalData::~CalData()
{

}

CalData::CalData(const CalData &r)
{
	goodness=r.goodness;
	IsoNum=r.IsoNum;
	mzt=r.mzt;
	scannum=r.scannum;
	tic=r.tic;
	ch=r.ch;
	RT=r.RT;
	//addtional 2009.10.10
	BaseInt=r.BaseInt;
	int i;
	for(i=0;i<IsoNum;i++)
	{
		IsotopicE[i]=r.IsotopicE[i];
		IsotopicT[i]=r.IsotopicT[i];
		Isomz[i]=r.Isomz[i];
		S_N[i]=r.S_N[i];
	}
	//防止意外赋值
	for(;i<MAX_ISO;i++)
	{
		IsotopicE[i]=0;
		IsotopicT[i]=0;
		Isomz[i]=0;
		S_N[i]=0;
	}
	//addtional on 2010.10.11
	size_t parnum=r.FTStatus.size();
	for(i=0;i<parnum;i++)	FTStatus.push_back(r.FTStatus[i]);
	//IIt=r.IIt;
	//ESt=r.ESt;
	//InstrumentType=r.InstrumentType;//default is unknown
}

CalData &CalData::operator=(const CalData &r)
{
	goodness=r.goodness;
	IsoNum=r.IsoNum;
	mzt=r.mzt;
	scannum=r.scannum;
	tic=r.tic;
	ch=r.ch;
	RT=r.RT;
	//addtional 2009.10.10
	BaseInt=r.BaseInt;
	int i;
	for(i=0;i<IsoNum;i++)
	{
		IsotopicE[i]=r.IsotopicE[i];
		IsotopicT[i]=r.IsotopicT[i];
		Isomz[i]=r.Isomz[i];
		S_N[i]=r.S_N[i];
	}
	//防止意外赋值
	for(;i<MAX_ISO;i++)
	{
		IsotopicE[i]=0;
		IsotopicT[i]=0;
		Isomz[i]=0;
		S_N[i]=0;
	}
	//addtional on 2010.10.11
	FTStatus.clear();
	size_t parnum=r.FTStatus.size();
	for(i=0;i<parnum;i++)	FTStatus.push_back(r.FTStatus[i]);
	//IIt=r.IIt;
	//ESt=r.ESt;
	//InstrumentType=r.InstrumentType;//default is unknown
	return *this;
}

double CalData::GetmaxInt() const
{
	int i;
	double maxint=0;
	for(i=0;i<MAX_ISO;i++)
	{
		if(IsotopicE[i]<=1e-2) break;
		if(maxint<IsotopicE[i]) maxint=IsotopicE[i];
	}
	return maxint;
}

double CalData::GetmaxIntR() const
{
	int i;
	double maxint=0;
	if(BaseInt<1e-2) return 0;
	for(i=0;i<MAX_ISO;i++)
	{
		if(IsotopicE[i]<=1e-2) break;	
		if(IsotopicE[i]>maxint) maxint=IsotopicE[i];
	}
	return maxint/BaseInt;
}
//this method is not fit for overlap peaks
int CalData::GetIsoNum()
{
	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		if(IsotopicE[i]<1e-2) break;
	}
	return i;
}

//fpr overlap peaks, can not calculate goodness using this method
double CalData::IsoMatchGD()
{
	int i,numMatch=0;	
	for(i=0;i<MAX_ISO;i++)
	{
		if(IsotopicE[i]>1e-2)numMatch++;
	}	
	double Intsum=0;
	double Intsum1=0;
	double Intsum2=0;
	for(i=0;i<MAX_ISO;i++)
	{
		Intsum+=IsotopicE[i]*IsotopicE[i];
		Intsum1+=IsotopicT[i]*IsotopicT[i];
		Intsum2+=IsotopicT[i]*IsotopicE[i];
	}
	Intsum=sqrt(Intsum1*Intsum);
	if(Intsum<=1e-6) return 0;	
	else return Intsum2/Intsum;
}

void CalData::OutPut2Str(CString &str,vector<size_t> &ParSel)
{
	int iso_num=GetIsoNum();
	double ty=(Isomz[0]-mzt)*1e6/Isomz[0];
	for(int i=0;i<iso_num;i++)
	{
		str.Format("%lf\t1:%lf\t2:%d\t3:%lf\t4:%lf\t5:%d\t6:%lf\t7:%lf",ty,Isomz[i],iso_num,RT,S_N[i],i,IsotopicE[i],tic);
		CString item;
		size_t parnum=ParSel.size();
		for(int i=0;i<parnum;i++)
		{
			item.Format("\t%d:%lf",i+8,FTStatus[ParSel[i]]);
			str+=item;
		}
		str+="\n";
	}
}


bool CalData::PackagePar(double *fea,int iso,vector<size_t> &ParSel) //package the SVM input vector
{	
	fea[0]=Isomz[iso];
	fea[1]=RT;
	fea[2]=tic;
	fea[3]=log10(IsotopicE[iso]+1e-4);
	fea[4]=sqrt(S_N[iso]);
	fea[5]=IsotopicE[iso]/BaseInt;
	size_t sg=ParSel.size();
	for(size_t i=0;i<sg;i++)
	{
		fea[FIXED_PAR_NUM_LIN+i]=FTStatus[ParSel[i]];
	}
	return true;
}

void CalData::PackagePar(double *fea,vector<size_t> &ParSel) //package the SVM input vector
{
	fea[0]=Isomz[0];
	fea[1]=RT;
	fea[2]=tic;
	fea[3]=log10(IsotopicE[0]+1e-4);
	fea[4]=sqrt(S_N[0]);
	fea[5]=IsotopicE[0]/BaseInt;
	size_t sg=ParSel.size();
	for(size_t i=0;i<sg;i++)
	{
		fea[FIXED_PAR_NUM_LIN+i]=FTStatus[ParSel[i]];
	}	
}

void CalData::PackageParFixed(double *fea,vector<size_t> &ParSel) //package the SVM input vector
{	
	fea[1]=RT;
	fea[2]=tic;
	size_t sg=ParSel.size();
	for(size_t i=0;i<sg;i++)
	{
		fea[FIXED_PAR_NUM_LIN+i]=FTStatus[ParSel[i]];
	}	
}

bool CalData::PackageParVal(double *fea,int iso) //package the SVM input vector
{
	fea[0]=Isomz[iso];
	fea[3]=log10(IsotopicE[iso]+1e-4);
	fea[4]=sqrt(S_N[iso]);
	fea[5]=IsotopicE[iso]/BaseInt;
	return true;
}


bool CalData::PackagePar(svm_node *fea,int iso,vector<size_t> &ParSel) //package the SVM input vector
{	
	fea[0].value=Isomz[iso]-MASS_P*iso/ch;
	fea[1].value=RT;
	fea[2].value=tic;
	fea[3].value=log10(IsotopicE[iso]+1e-4);
	fea[4].value=sqrt(S_N[iso]);
	fea[5].value=IsotopicE[iso]/BaseInt;	
	size_t sg=ParSel.size();
	for(size_t i=0;i<sg;i++)
	{
		fea[FIXED_PAR_NUM_SVM+i].value=FTStatus[ParSel[i]];
	}
	return true;
}

void CalData::PackagePar(svm_node *fea,vector<size_t> &ParSel) //package the SVM input vector
{
	fea[0].value=Isomz[0];
	fea[1].value=RT;
	fea[2].value=tic;
	fea[3].value=log10(IsotopicE[0]+1e-4);
	fea[4].value=sqrt(S_N[0]);	
	fea[5].value=IsotopicE[0]/BaseInt;	
	size_t sg=ParSel.size();
	for(size_t i=0;i<sg;i++)
	{
		fea[FIXED_PAR_NUM_SVM+i].value=FTStatus[ParSel[i]];
	}	
}

void CalData::PackageParFixed(svm_node *fea,vector<size_t> &ParSel) //package the SVM input vector
{	
	fea[1].value=RT;
	fea[2].value=tic;
	size_t sg=ParSel.size();
	for(size_t i=0;i<sg;i++)
	{
		fea[FIXED_PAR_NUM_SVM+i].value=FTStatus[ParSel[i]];
	}	
}

bool CalData::PackageParVal(svm_node *fea,int iso) //package the SVM input vector
{
	fea[0].value=Isomz[iso]-MASS_P*iso/ch;
	fea[3].value=log10(IsotopicE[iso]+1e-4);
	fea[4].value=sqrt(S_N[iso]);
	fea[5].value=IsotopicE[iso]/BaseInt;	;
	return true;
}



bool CalData::IspmzIn(double pmz,double dmppm)
{
	int i;
	for(i=0;i<MAX_ISO;i++)
	{
		if(IsotopicE[i]<=1e-2) break;
		double dm=pmz-Isomz[i];
		if(fabs(dm)<pmz*dmppm/1e6) return true;
	}
	return false;
}

//end of define CalData
//////////////////////////////////////////////////////////////////////
CReport::CReport()
{	
	count=0;
	mcount=0;
	IsRepValidate=false;
	Error_code=0x00;
	for(int i=0;i<4;i++)METModel[i]=0;
	IsMETModel=false;
	Model_mean=0;
	Model_std=0;
	Model_gd=0;
	ppbLevel=0;
}

CReport::~CReport()
{

}

CReport::CReport(const CReport &r)
{	
	rFile=r.rFile;
	oPath=r.oPath;
	RawFile=r.RawFile;	
	count=r.count;
	mcount=r.mcount;
	mtA=r.mtA;
	mtB=r.mtB;
	IsRepValidate=r.IsRepValidate;
	Error_code=r.Error_code;
	for(int i=0;i<4;i++)METModel[i]=r.METModel[i];
	IsMETModel=r.IsMETModel;
	Model_mean=r.Model_mean;
	Model_std=r.Model_std;
	Model_gd=r.Model_gd;
	size_t i, sg=r.SelParName.size();
	for(i=0;i<sg;i++)
	{
		SelParName.push_back(r.SelParName[i]);
	}

	sg=r.LCoff.size();
	for(i=0;i<sg;i++)
	{
		LCoff.push_back(r.LCoff[i]);
	}	

	ppbLevel=r.ppbLevel;
}


CReport &CReport::operator=(const CReport &r)
{	
	rFile=r.rFile;
	oPath=r.oPath;
	RawFile=r.RawFile;	
	count=r.count;
	mcount=r.mcount;
	mtA=r.mtA;
	mtB=r.mtB;
	IsRepValidate=r.IsRepValidate;
	Error_code=r.Error_code;
	for(int i=0;i<4;i++)METModel[i]=r.METModel[i];
	IsMETModel=r.IsMETModel;
	Model_mean=r.Model_mean;
	Model_std=r.Model_std;
	Model_gd=r.Model_gd;
	SelParName.clear();
	size_t i,sg=r.SelParName.size();
	for(i=0;i<sg;i++)
	{
		SelParName.push_back(r.SelParName[i]);
	}

	LCoff.clear();
	sg=r.LCoff.size();
	for(i=0;i<sg;i++)
	{
		LCoff.push_back(r.LCoff[i]);
	}	

	ppbLevel=r.ppbLevel;
	return *this;
}

///common used functions
void GetIsoDis(double mass,double IsoDisT[MAX_ISO])
{	
	if(mass<=0.1)
	{
		IsoDisT[0]=1.0;
		return;
	}
	double tIsoDisT[6];

	tIsoDisT[0]=1.007*exp(-0.0005792*mass);

	double item=0.0006321*mass-0.09212 ;
	tIsoDisT[1]=item*exp(-item);

	item=0.0005683 *mass+0.02292;
	tIsoDisT[2]=item*item*exp(-item)/2;

	item=0.0005526*mass+0.09675;
	tIsoDisT[3]=item*item*item*exp(-item)/6;
	
	item=0.000568*mass+0.1138;
	tIsoDisT[4]=item*item*item*item*exp(-item)/24;
	
	item=0.0005795*mass+0.1215;
	tIsoDisT[5]=item*item*item*item*item*exp(-item)/120;
	double sum=0;
	int i;
	if(MAX_ISO<=6)
	{
		for(i=0;i<MAX_ISO;i++) 
		{
			sum+=tIsoDisT[i];
			IsoDisT[i]=tIsoDisT[i];
		}
		for(i=0;i<MAX_ISO;i++) IsoDisT[i]/=sum;
	}
	else
	{
		for(i=0;i<6;i++) 
		{
			sum+=tIsoDisT[i];
			IsoDisT[i]=tIsoDisT[i];
		}
		double tit=item*item*item*item*item;
		int ct=120;
		for(i=6;i<MAX_ISO;i++)
		{
			tit*=item;
			ct*=i;
			IsoDisT[i]=tit*exp(-item)/ct;
			sum+=IsoDisT[i];
		}		
		for(i=0;i<MAX_ISO;i++) IsoDisT[i]/=sum;
	}
}

void GetIsoDis(double mass,double *IsoDisT,int max_iso)
{	
	if(mass<=0.1)
	{
		IsoDisT[0]=1.0;
		return;
	}
		
	double tIsoDisT[6];

	tIsoDisT[0]=1.007*exp(-0.0005792*mass);

	double item=0.0006321*mass-0.09212 ;
	tIsoDisT[1]=item*exp(-item);

	item=0.0005683 *mass+0.02292;
	tIsoDisT[2]=item*item*exp(-item)/2;

	item=0.0005526*mass+0.09675;
	tIsoDisT[3]=item*item*item*exp(-item)/6;
	
	item=0.000568*mass+0.1138;
	tIsoDisT[4]=item*item*item*item*exp(-item)/24;
	
	item=0.0005795*mass+0.1215;
	tIsoDisT[5]=item*item*item*item*item*exp(-item)/120;
	double sum=0;
	int i;
	if(max_iso<=6)
	{
		for(i=0;i<max_iso;i++) 
		{
			sum+=tIsoDisT[i];
			IsoDisT[i]=tIsoDisT[i];
		}
		for(i=0;i<max_iso;i++) IsoDisT[i]/=sum;
	}
	else
	{
		for(i=0;i<6;i++) 
		{
			sum+=tIsoDisT[i];
			IsoDisT[i]=tIsoDisT[i];
		}
		double tit=item*item*item*item*item;
		int ct=120;
		for(i=6;i<max_iso;i++)
		{
			tit*=item;
			ct*=i;
			IsoDisT[i]=tit*exp(-item)/ct;
			sum+=IsoDisT[i];
		}		
		for(i=0;i<max_iso;i++) IsoDisT[i]/=sum;
	}
}

double GetIsoDis(double mass,int iso_no)
{
	double item;
	switch(iso_no)
	{
	case 1:
		return 1.007*exp(-0.0005792*mass);
		break;
	case 2:
		item=0.0006321*mass-0.09212 ;
		return item*exp(-item);
		break;
	case 3:
		item=0.0005683 *mass+0.02292;
		return item*item*exp(-item)/2;
		break;
	case 4:
		item=0.0005526*mass+0.09675;
		return item*item*item*exp(-item)/6;
		break;
	case 5:
		item=0.000568*mass+0.1138;
		return item*item*item*item*exp(-item)/24;
		break;
	case 6:
		item=0.0005795*mass+0.1215;
		return item*item*item*item*item*exp(-item)/120;
		break;
	default:
		item=0.0005795*mass+0.1215;
		double tit=item*item*item*item*item;
		int ct=120;
		for(int i=6;i<=iso_no;i++)
		{
			tit*=item;
			ct*=i;	
		}	
		return tit*exp(-item)/ct;
	}
}


bool bootstrp(double *x,int num,double *output,int numout,statsF pF)
{
	if(num<=0||numout<=0) return false;
	double *ty;
	ty=new double[num];
	srand((unsigned)time(NULL));
	for(int i=0;i<numout;i++)
	{
		for(int k=0;k<num;k++)
		{		
			int r=(int)(rand()*1.0*num/RAND_MAX);
			if(r>=num)r=num-1;
			ty[k]=x[r];
		}
		output[i]=pF(ty,num);
	}
	delete []ty;
	return true;
}

bool bootstrp(vector<double> &x,vector<double> &y,int numout,statsF pF)
{
	size_t num=x.size();
	y.clear();
	if(num<=0||numout<=0) return false;
	double *ty;
	ty=new double[num];
	srand((unsigned)time(NULL));
	for(int i=0;i<numout;i++)
	{
		for(int k=0;k<num;k++)
		{		
			int r=(int)(rand()*1.0*num/RAND_MAX);
			if(r>=num)r=num-1;
			ty[k]=x[r];
		}
		y.push_back(pF(ty,num));
	}
	delete []ty;
	return true;

}

Pre_calData::Pre_calData()
{
	mzt=0;
	mze=0;
	charge=0;
	scan=0;

}

Pre_calData::~Pre_calData()
{

}

Pre_calData::Pre_calData(const Pre_calData &cp)
{
	mzt=cp.mzt;
	mze=cp.mze;
	charge=cp.charge;
	scan=cp.scan;
}

Pre_calData& Pre_calData::operator=(const Pre_calData &cp)
{
	mzt=cp.mzt;
	mze=cp.mze;
	charge=cp.charge;
	scan=cp.scan;
	return *this;
}

//default advance parameters
AdvPars::AdvPars()
{
	Min_ISO=1;
	Min_Rel_Int=0.001;
	Iso_GD=0.2;
	Xic_Int_Max=3;
	RT_Max_Int=0.2;
	Ch_RT_Min=-5;
	Ch_RT_Max=5;
	PIonResearch=true;
	XICCalibrate=true;
	TrainSize=6000;
	Pre_dm=0;//ppm,means use auto
	SVM_C=64;
	Max_Ch=4;
	MinIsoNum=2;
	RIntCTH=0.01;
	GdCtH=0.2;
	LExt=10;
	RExt=10;
	//add new
	MinRT=-5;
	MaxRT=5;
	//add new
	Min_SN=2;

	XIC_SW=21;
	XIC_SMC=7;
	IsDyaXIC=false;
}

AdvPars::~AdvPars()
{

}

AdvPars::AdvPars(const AdvPars &cp)
{
	Min_ISO=cp.Min_ISO;
	Min_Rel_Int=cp.Min_Rel_Int;
	Iso_GD=cp.Iso_GD;
	Xic_Int_Max=cp.Xic_Int_Max;
	RT_Max_Int=cp.RT_Max_Int;
	Ch_RT_Min=cp.Ch_RT_Min;
	Ch_RT_Max=cp.Ch_RT_Max;
	PIonResearch=cp.PIonResearch;
	XICCalibrate=cp.XICCalibrate;
	TrainSize=cp.TrainSize;
	Pre_dm=cp.Pre_dm;//ppm	
	SVM_C=cp.SVM_C;

	Max_Ch=cp.Max_Ch;
	MinIsoNum=cp.MinIsoNum;
	RIntCTH=cp.RIntCTH;
	GdCtH=cp.GdCtH;
	LExt=cp.LExt;
	RExt=cp.RExt;
	MinRT=cp.MinRT;
	MaxRT=cp.MaxRT;

	Min_SN=cp.Min_SN;

	XIC_SW=cp.XIC_SW;
	XIC_SMC=cp.XIC_SMC;
	IsDyaXIC=cp.IsDyaXIC;
}

AdvPars &AdvPars::operator=(const AdvPars &cp)
{
	Min_ISO=cp.Min_ISO;
	Min_Rel_Int=cp.Min_Rel_Int;
	Iso_GD=cp.Iso_GD;
	Xic_Int_Max=cp.Xic_Int_Max;
	RT_Max_Int=cp.RT_Max_Int;
	Ch_RT_Min=cp.Ch_RT_Min;
	Ch_RT_Max=cp.Ch_RT_Max;
	PIonResearch=cp.PIonResearch;
	XICCalibrate=cp.XICCalibrate;
	TrainSize=cp.TrainSize;
	Pre_dm=cp.Pre_dm;//ppm	
	SVM_C=cp.SVM_C;

	Max_Ch=cp.Max_Ch;
	MinIsoNum=cp.MinIsoNum;
	RIntCTH=cp.RIntCTH;
	GdCtH=cp.GdCtH;
	LExt=cp.LExt;
	RExt=cp.RExt;

	MinRT=cp.MinRT;
	MaxRT=cp.MaxRT;

	Min_SN=cp.Min_SN;

	XIC_SW=cp.XIC_SW;
	XIC_SMC=cp.XIC_SMC;
	IsDyaXIC=cp.IsDyaXIC;
	return *this;
}

bool AdvPars::CheckValidate()
{
	if(Min_ISO<1||Min_ISO>10) return false;
	if(Min_Rel_Int<0||Min_Rel_Int>0.99999) return false;
	if(Iso_GD<0||Iso_GD>0.99999) return false;
	if(Xic_Int_Max<0||Xic_Int_Max>100) return false;
	if(RT_Max_Int<=0||RT_Max_Int>20) return false;
	if(Ch_RT_Min>=0||Ch_RT_Min<=-20) return false;
	if(Ch_RT_Max<=0||Ch_RT_Max>=20) return false;
	if(TrainSize<=300) return false;
	//if(Pre_dm<=1||Pre_dm>100) return false;
	if(SVM_C<=0||SVM_C>1024) return false;
	if(Max_Ch<=0||Max_Ch>20) return false;
	if(MinIsoNum<=1||MinIsoNum>10) return false;
	if(RIntCTH<=1e-4||RIntCTH>=0.99999) return false;
	if(GdCtH<=0.1||GdCtH>=0.99999) return false;
	if(LExt<0||LExt>1000) return false;
	if(RExt<0||RExt>1000) return false;
	if(MinRT>0||MinRT<-1000) return false;
	if(MaxRT<0||MaxRT>1000) return false;
	//if(Min_SN<0） return false;
	if(XIC_SW<7) return false;
	if(XIC_SMC<3) return false;	
	return true;
}

void AdvPars::SaveToFile(CStdioFile &dFile)
{
	CString tmp;
	tmp.Format("Min_ISO=%d\n",Min_ISO);
	dFile.WriteString(tmp);

	tmp.Format("Min_Rel_Int=%lf\n",Min_Rel_Int);
	dFile.WriteString(tmp);

	tmp.Format("Iso_GD=%lf\n",Iso_GD);
	dFile.WriteString(tmp);

	tmp.Format("Xic_Int_Max=%d\n",Xic_Int_Max);
	dFile.WriteString(tmp);

	tmp.Format("RT_Max_Int=%lf\n",RT_Max_Int);
	dFile.WriteString(tmp);

	tmp.Format("Ch_RT_Min=%lf\n",Ch_RT_Min);
	dFile.WriteString(tmp);
	
	tmp.Format("Ch_RT_Max=%lf\n",Ch_RT_Max);
	dFile.WriteString(tmp);

	tmp.Format("PIonResearch=%d\n",PIonResearch);
	dFile.WriteString(tmp);

	tmp.Format("XICCalibrate=%d\n",XICCalibrate);
	dFile.WriteString(tmp);

	tmp.Format("TrainSize=%d\n",TrainSize);
	dFile.WriteString(tmp);
	
	tmp.Format("Pre_dm=%lf\n",Pre_dm);
	dFile.WriteString(tmp);

	tmp.Format("SVM_C=%lf\n",SVM_C);
	dFile.WriteString(tmp);

	tmp.Format("Max_Ch=%d\n",Max_Ch);
	dFile.WriteString(tmp);

	tmp.Format("MinIsoNum=%d\n",MinIsoNum);
	dFile.WriteString(tmp);

	tmp.Format("RIntCTH=%lf\n",RIntCTH);
	dFile.WriteString(tmp);

	tmp.Format("GdCtH=%lf\n",GdCtH);
	dFile.WriteString(tmp);

	tmp.Format("LExt=%d\n",LExt);
	dFile.WriteString(tmp);

	tmp.Format("RExt=%d\n",RExt);
	dFile.WriteString(tmp);

	tmp.Format("RT_Min=%lf\n",MinRT);
	dFile.WriteString(tmp);
	
	tmp.Format("RT_Max=%lf\n",MaxRT);
	dFile.WriteString(tmp);

	tmp.Format("Min_SN=%lf\n",Min_SN);
	dFile.WriteString(tmp);

	tmp.Format("SW=%d\n",XIC_SW);
	dFile.WriteString(tmp);

	tmp.Format("SMC=%d\n",XIC_SMC);
	dFile.WriteString(tmp);

	tmp.Format("IsDyaXIC=%d\n",IsDyaXIC);
	dFile.WriteString(tmp);

}

bool AdvPars::LoadFromFile(CStdioFile &dFile)
{
	CString tmp;

	if(!dFile.ReadString(tmp)) return false;
	int RD=sscanf((LPCTSTR)tmp,"Min_ISO=%d",&Min_ISO);
	if(RD!=1) return false;	

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Min_Rel_Int=%lf",&Min_Rel_Int);
	if(RD!=1) return false;		

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Iso_GD=%lf",&Iso_GD);
	if(RD!=1) return false;	

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Xic_Int_Max=%d",&Xic_Int_Max);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"RT_Max_Int=%lf",&RT_Max_Int);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Ch_RT_Min=%lf",&Ch_RT_Min);
	if(RD!=1) return false;
	
	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Ch_RT_Max=%lf",&Ch_RT_Max);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"PIonResearch=%d",&PIonResearch);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"XICCalibrate=%d",&XICCalibrate);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"TrainSize=%d",&TrainSize);
	if(RD!=1) return false;
	
	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Pre_dm=%lf",&Pre_dm);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"SVM_C=%lf",&SVM_C);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Max_Ch=%d",&Max_Ch);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"MinIsoNum=%d",&MinIsoNum);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"RIntCTH=%lf",&RIntCTH);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"GdCtH=%lf",&GdCtH);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"LExt=%d",&LExt);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"RExt=%d",&RExt);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"RT_Min=%lf\n",&MinRT);
	if(RD!=1) return false;
		
	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"RT_Max=%lf\n",&MaxRT);	
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"Min_SN=%lf\n",&Min_SN);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"SW=%d\n",&XIC_SW);
	if(RD!=1) return false;
	
	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"SMC=%d\n",&XIC_SMC);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"IsDyaXIC=%d\n",&IsDyaXIC);
	if(RD!=1) return false;

	return true;
}