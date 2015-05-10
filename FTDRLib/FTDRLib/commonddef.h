#pragma once
#include "afxmt.h"
#include "MResult.h"
#include "OutRecord.h"
#include "HistList.h"
#include <string>
#include <vector>
#include "MSTypes.h"
#include "svm.h"
using std::vector;
using std::string;

extern int ConselRow;//保存当前调试窗口光标的行位置
extern int ConselCol;//保存单签调试窗口光标的列位置

//post message
#define MSG_THREAD_END NULL
#define MSG_TASK_FINISHED 0x01
#define MSG_RAW_LOAD 0x02
#define MSG_RESULT_LOAD 0x04
#define MSG_MODEL_END 0x08
#define MSG_NOIDRESULTS 0x10

//error message code
#define ERROR_RAWFAILUR 0x01//open raw file failure
#define ERROR_NOTLOAD 0x02//no records for model building loaded
#define ERROR_COUPUTF 0x04 //create output file failure
#define ERROR_MDFAILUR 0x08//calibrating model building failure, too less validated identifications for model building
#define ERROR_RNOPEN 0x10//can notopen the results file or path no out files
#define ERROR_NOS_RFORMAT 0x20//the extend calibration model require the thermo raw input, if not return this message
#define ERROR_MS1OUTF 0x40
#define ERROR_MODEL 0x80


//perform isotopic calibration
#define MAX_ISO 6//最多考虑的同位素峰数量

//define the result type
#define R_TYPE_MASCOT 0//Mascot
#define R_TYPE_SEQUEST 1//Sequest
#define R_TYPE_PEPXML 2//PepXML
#define R_TYPE_MASSMATRIX 3//massmatrix

#define TIC_SCALE 1e7f//TIC归一化因子
#define MZE_SCALE 1e3f//质荷比归一化因子
#define PINT_SCALE 1e4//母离子信号强度归一化因子

extern double MINRT;//普通XIC搜索中，RT范围最小值，默认-5.0min
extern double MAXRT;//普通XIC搜索中，RT范围最大值，默认5.0min
extern double GD_CUT;//=0.6;同位素峰的你和优度门限
extern double MIN_SN;//信噪比门限，默认为2
extern int MAX_CH;//=4;最大的考虑电荷
//extern int IT_LIM;//=5; not used
extern double RINT_CUT;//=0.001;相对信号强度过滤门限
#define MAX_TITLE 256//title的字符串长度，在mgf文件写入时使用
//#define MASS_LIM 300//not used
extern int ISO_CUT;//=1;//means no cut here同位素峰数目的限制
extern int INT_TIME_MAX;// 3//XIC构建中，信号缺失的最大次数
extern double MAX_INT_RT;// 0.2//XIC截断中，允许的时间最大间隔
#define CH_UNKNOWN 0//未知电荷

extern double CH_RT_MIN;//=-5;电荷确定中，搜索图谱叠加的区间向左扩展范围
extern double CH_RT_MAX;//=5;电荷确定中，搜索图谱叠加的区间向右扩展范围

extern long MAX_SVM_TRAIN;// 10000 SVM模型训练样本数量最大值

extern int MIN_H_ISO_NUM; //=2扩展查找母离子时，同位素峰最小值
extern double R_INT_CUT_H;//= 0.01扩展查找时，相对信号强度最小值
extern double GD_CUT_H;//= 0.2 扩展查找时，同位素峰拟合优度门限

#define MASS_H 1.007825 //mass for Hydrogen
#define MASS_H_ME 1.007276 //mass for Hydrogen without Electron
#define MASS_E 0.000549 //mass for Electron
#define MASS_P 1.0033548

extern int LLEXTEND;// 10
extern int RLEXTEND;// 10

extern double PRE_DM;//10ppm
extern double PRE_DM_10;
extern bool XICCALIBRATE;//true
extern bool PARIONRED;//true;
extern double SVMC;

extern int DYA_XIC_SW;
extern int DYA_XIC_SM;
extern bool IS_DYA_XIC;

#define GD_CUT_L 0.01//在group拟合函数中使用的过滤，为了保证原来的母离子不丢失，这个过滤门限仅仅选择为0.1
#define PPMLIM_MIN -30
#define PPMLIM_MAX 30

#define MODEL_LIN 0x01
#define MODEL_SVM 0x02
#define MODEL_NON 0x00
#define MODEL_ERROR -1

#define MS2NON 0x01
#define MS2MGF 0x02
#define MS2EMGF 0x03

#define SQRTTWO 1.4142135623730950488016887242097

//define the isotopic mass difference
const static double ISO_DIFF[MAX_ISO]={
	/*0.0,统计参数
	1.003308559119778,
	2.0062822183398778,
	3.007692024454673,
	4.009119148582318,
	5.010832505631107*/
	//理论参数
	0.0,
	1.0033548,//mass of neutron
	2.0067096,
	3.0100644,
	4.0134192,
	5.0167740
};
//to support the status and finial report
class CReport
{
public:
	CReport();
	~CReport();
	CString rFile;//搜库结果文件完整路径
	CString RawFile;//raw文件完整路径
	CString oPath;//输出路径
	int count;//建模使用的结果数量
	int mcount;//校正结果数量
	bool IsRepValidate;//内部使用，确定该报告是否有效
	myHistogram mtA;//用于显示，校正前母离子误差分布直方图
	myHistogram mtB;//用于显示，校正后母离子误差分布直方图
	double METModel[4];//信号强度相关误差分布模型参数
	bool IsMETModel;//是否使用信号强度相关模型，确定上述参数是否有效
	int Error_code;//错误代码，已经预定义
	CReport(const CReport &r);
	CReport &operator=(const CReport &r);
public:
	double Model_mean;//svm或者Lin训练残差均值
	double Model_std;//svm或者Lin训练残差标准差
	double Model_gd;//SVM或者LIn模型拟合优度：RSquare
public:
	vector<string> SelParName;//特征选择结果
	vector<double> LCoff;//线性模型的系数, 和上述参数对应
public:
	int ppbLevel;//ppb级别数量统计
};

#define MS1_NON 0x00
#define MS1_mzML 0x01//输出一级图谱类型
#define MS1_mzXML 0x02//
//trans the parameters between the GUI thread and work threads
class ThreadParm
{
private://维持的报告列表，由线程数决定
	CReport *RepList;
	int count;
public://传递参数
	//参数读取锁定信号
	CCriticalSection GL_CriticalSection;
	//同步事件
	HANDLE hEventP;//任务是否准备好的状态指示，在workThread中使用
	HANDLE hEventE;//任务参数是不是已经被使用的状态指示，在workThread中设置，在CFTDRDlg类中判断
	HANDLE hEventB;//没有使用
	//搜库结果文件
	CString RFile;
	//Raw文件
	CString RawFile;
	//输出路径，输出校正结果，模型和报告
	CString sOutPath;
	//SVM模型文件
	//CString SVMModelFile;
	//线程ID
	UINT ThreadID;
	//当前任务是不是结束，同步使用
	bool IsEnd;
	//搜库结果过滤参数，Sequest和mascot分别处理
	//如果使用类继承方式，还可以实现多态，只需要一个基类的指针就可以了
	//这里是可以改进的地方
	Filter dFilter;
	SFilter dFilter1;
	//搜库结果类型
	int type;	
	//该任务是不是已经处理过
	bool IsDeal;
	//是否输出SVM模型
	bool IsOPSVMModel;
	//是否输出校正模型建模数据
	bool IsOpCData;
	//是不是输出MS1数据，输出格式是什么
	int ms1format;
	//MS2输出格式，包括不输出
	int ms2format;
	//是否使用FT状态参数
	bool IsUseStatsPar;
	//校正模型选择，包括不用任何模型，线性模型和svm模型3种
	int ModelType;
public:
	ThreadParm();
	ThreadParm(const ThreadParm &r);
	~ThreadParm();	
	//分配/释放报告产生的存储空间
	void AllocRepBuff(int n);
	void FreeRepBuff();
	//填写/读取报告内容
	bool SetRepBuff(int idx,CReport *rt);
	bool GetRepBuff(int idx,CReport *rt);
	//重置报告状态
	void ResetRepState(int idx);
	ThreadParm &operator=(const ThreadParm &r);
};

//高级参数设置
class AdvPars
{
public:
	int Min_ISO;//最小同位素峰数目
	double Min_Rel_Int;//最小相对信号强度
	double Iso_GD;//最小拟合优度
	int Xic_Int_Max;//XIC搜索中最大gap次数，超过该次数就截断
	double RT_Max_Int;//XIC搜索允许最大时间间隔
	double Ch_RT_Min;//电荷确定的区间最小值
	double Ch_RT_Max;//电荷确定的区间最小值
	double Pre_dm;//预定义母离子搜索误差，在XIC搜索中使用
	//double IsoDiff[MAX_ISO];
	//double stats_Iso_MET[MAX_ISO];
	double SVM_C;//svm训练参数，越大需要时间越长
	//double SVM_gama;
	bool PIonResearch;//是否进行母离子扩展搜索
	bool XICCalibrate;//是否执行XIC校正
	long TrainSize;//svm训练样本最大数量，如果给定样本大于该数值，进行随机选择
	int Max_Ch;//最大的考虑电荷
	int MinIsoNum;//扩展查找母离子时，同位素峰最小值
	double RIntCTH;//扩展查找时，相对信号强度最小值
	double GdCtH;
	int LExt;
	int RExt;
	//to limit the search work
	double MinRT;
	double MaxRT;
	//
	double Min_SN;//最低信噪比
	//动态XIC截断参数
	int XIC_SW;
	int XIC_SMC;
	bool IsDyaXIC;
	//
	AdvPars();
	~AdvPars();
	AdvPars(const AdvPars &cp);
	AdvPars &operator=(const AdvPars &cp);
	bool CheckValidate();
	void SaveToFile(CStdioFile &dFile);
	bool LoadFromFile(CStdioFile &dFile);
};

#define FT_PAR_NUM 7
#define Repeat_IDX 5//there are two parameters with same name, we slected the 2nd one
const static string FT_Status_Par[FT_PAR_NUM]=
{
	"Source Current (uA)",
	"RF Detector Temp (C)",
	"RF Generator Temp (C)",	
	"Ambient Temp. (C)",
	"FT EA Temp. (C)",
	"FT RF1 Amp. Temp. (C)",	
	"Nitrogen (%)"
};

//#define INS_FT 0x01
//#define INS_ORBITRAP 0x02
#define MAX_FEA_NUM 12
#define FIXED_PAR_NUM_SVM 6
#define FIXED_PAR_NUM_LIN 6
//one associate with a isotopic cluster, because they share too many parameters
class CalData
{
public:
	//basic parameters for the data loaded
	double goodness;
	int IsoNum;
	double tic;
	double mzt;
	int ch;
	long scannum;
	double RT;
	//addtional parameters, added on 2009.10.9
	double BaseInt;
	double IsotopicE[MAX_ISO];
	double IsotopicT[MAX_ISO];
	double Isomz[MAX_ISO];
	//added on 2011.6.2
	double S_N[MAX_ISO];
	//End of this section
	//****************the following modified on 2010.10.11 by zhangjy****************************//
	//addtional parameters for FT status in raw file,added on 2010.10.11
	//read by function GetStatusLogForScanNum 
	vector<double> FTStatus;//modified on 2010.12.3, if in calibration, not all is used
	//addtional parameters for Ion Injection Time and Elapsed Scan Time in raw file,added on 2010.10.11
	//read by function GetTrailerExtraForScanNum
	//double IIt;//Orbitrap needed
	//double ESt;//Orbitrap needed //merge to the FTStatus on 2011.11.8
	//MSInstrumentModelType InstrumentType;
public:
	CalData();
	~CalData();
	CalData(const CalData &r);
	int GetIsoNum();
	CalData &operator=(const CalData &r);
	void OutPut2Str(CString &str,vector<size_t> &ParSel);//output all the calibrate features
	double IsoMatchGD();//get the isotopic map match goodness
	//for generally use
	bool PackagePar(double *fea,int iso,vector<size_t> &ParSel); //package the SVM input vector	
	void PackageParFixed(double *fea,vector<size_t> &ParSel); //package the SVM input vector
    bool PackageParVal(double *fea,int iso);//package the SVM input vector
	void PackagePar(double *fea,vector<size_t> &ParSel); //package the SVM input vector
	//for svm
	bool PackagePar(svm_node *fea,int iso,vector<size_t> &ParSel); //package the SVM input vector	
	void PackageParFixed(svm_node *fea,vector<size_t> &ParSel); //package the SVM input vector
    bool PackageParVal(svm_node *fea,int iso);//package the SVM input vector
	void PackagePar(svm_node *fea,vector<size_t> &ParSel); //package the SVM input vector		
	//
	bool IspmzIn(double pmz,double dmppm);
	double GetmaxInt() const;
	double GetmaxIntR() const;
};

class Pre_calData
{
public:
	double mzt;
	double mze;
	int charge;
	long scan;
public:
	Pre_calData();
	~Pre_calData();
	Pre_calData(const Pre_calData &cp);
	Pre_calData& operator=(const Pre_calData &cp);
};

//common functions
void GetIsoDis(double mass,double IsoDisT[MAX_ISO]);
void GetIsoDis(double mass,double *IsoDisT,int max_iso);
double GetIsoDis(double mass,int iso_no);

typedef double(*statsF)(double *x,int num);
bool bootstrp(double *x,int num,double *output,int numout,statsF pF);
bool bootstrp(vector<double> &x,vector<double> &y,int numout,statsF pF);
