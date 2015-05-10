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

extern int ConselRow;//���浱ǰ���Դ��ڹ�����λ��
extern int ConselCol;//���浥ǩ���Դ��ڹ�����λ��

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
#define MAX_ISO 6//��࿼�ǵ�ͬλ�ط�����

//define the result type
#define R_TYPE_MASCOT 0//Mascot
#define R_TYPE_SEQUEST 1//Sequest
#define R_TYPE_PEPXML 2//PepXML
#define R_TYPE_MASSMATRIX 3//massmatrix

#define TIC_SCALE 1e7f//TIC��һ������
#define MZE_SCALE 1e3f//�ʺɱȹ�һ������
#define PINT_SCALE 1e4//ĸ�����ź�ǿ�ȹ�һ������

extern double MINRT;//��ͨXIC�����У�RT��Χ��Сֵ��Ĭ��-5.0min
extern double MAXRT;//��ͨXIC�����У�RT��Χ���ֵ��Ĭ��5.0min
extern double GD_CUT;//=0.6;ͬλ�ط������Ŷ�����
extern double MIN_SN;//��������ޣ�Ĭ��Ϊ2
extern int MAX_CH;//=4;���Ŀ��ǵ��
//extern int IT_LIM;//=5; not used
extern double RINT_CUT;//=0.001;����ź�ǿ�ȹ�������
#define MAX_TITLE 256//title���ַ������ȣ���mgf�ļ�д��ʱʹ��
//#define MASS_LIM 300//not used
extern int ISO_CUT;//=1;//means no cut hereͬλ�ط���Ŀ������
extern int INT_TIME_MAX;// 3//XIC�����У��ź�ȱʧ��������
extern double MAX_INT_RT;// 0.2//XIC�ض��У������ʱ�������
#define CH_UNKNOWN 0//δ֪���

extern double CH_RT_MIN;//=-5;���ȷ���У�����ͼ�׵��ӵ�����������չ��Χ
extern double CH_RT_MAX;//=5;���ȷ���У�����ͼ�׵��ӵ�����������չ��Χ

extern long MAX_SVM_TRAIN;// 10000 SVMģ��ѵ�������������ֵ

extern int MIN_H_ISO_NUM; //=2��չ����ĸ����ʱ��ͬλ�ط���Сֵ
extern double R_INT_CUT_H;//= 0.01��չ����ʱ������ź�ǿ����Сֵ
extern double GD_CUT_H;//= 0.2 ��չ����ʱ��ͬλ�ط�����Ŷ�����

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

#define GD_CUT_L 0.01//��group��Ϻ�����ʹ�õĹ��ˣ�Ϊ�˱�֤ԭ����ĸ���Ӳ���ʧ������������޽���ѡ��Ϊ0.1
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
	/*0.0,ͳ�Ʋ���
	1.003308559119778,
	2.0062822183398778,
	3.007692024454673,
	4.009119148582318,
	5.010832505631107*/
	//���۲���
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
	CString rFile;//�ѿ����ļ�����·��
	CString RawFile;//raw�ļ�����·��
	CString oPath;//���·��
	int count;//��ģʹ�õĽ������
	int mcount;//У���������
	bool IsRepValidate;//�ڲ�ʹ�ã�ȷ���ñ����Ƿ���Ч
	myHistogram mtA;//������ʾ��У��ǰĸ�������ֲ�ֱ��ͼ
	myHistogram mtB;//������ʾ��У����ĸ�������ֲ�ֱ��ͼ
	double METModel[4];//�ź�ǿ��������ֲ�ģ�Ͳ���
	bool IsMETModel;//�Ƿ�ʹ���ź�ǿ�����ģ�ͣ�ȷ�����������Ƿ���Ч
	int Error_code;//������룬�Ѿ�Ԥ����
	CReport(const CReport &r);
	CReport &operator=(const CReport &r);
public:
	double Model_mean;//svm����Linѵ���в��ֵ
	double Model_std;//svm����Linѵ���в��׼��
	double Model_gd;//SVM����LInģ������Ŷȣ�RSquare
public:
	vector<string> SelParName;//����ѡ����
	vector<double> LCoff;//����ģ�͵�ϵ��, ������������Ӧ
public:
	int ppbLevel;//ppb��������ͳ��
};

#define MS1_NON 0x00
#define MS1_mzML 0x01//���һ��ͼ������
#define MS1_mzXML 0x02//
//trans the parameters between the GUI thread and work threads
class ThreadParm
{
private://ά�ֵı����б����߳�������
	CReport *RepList;
	int count;
public://���ݲ���
	//������ȡ�����ź�
	CCriticalSection GL_CriticalSection;
	//ͬ���¼�
	HANDLE hEventP;//�����Ƿ�׼���õ�״ָ̬ʾ����workThread��ʹ��
	HANDLE hEventE;//��������ǲ����Ѿ���ʹ�õ�״ָ̬ʾ����workThread�����ã���CFTDRDlg�����ж�
	HANDLE hEventB;//û��ʹ��
	//�ѿ����ļ�
	CString RFile;
	//Raw�ļ�
	CString RawFile;
	//���·�������У�������ģ�ͺͱ���
	CString sOutPath;
	//SVMģ���ļ�
	//CString SVMModelFile;
	//�߳�ID
	UINT ThreadID;
	//��ǰ�����ǲ��ǽ�����ͬ��ʹ��
	bool IsEnd;
	//�ѿ������˲�����Sequest��mascot�ֱ���
	//���ʹ����̳з�ʽ��������ʵ�ֶ�̬��ֻ��Ҫһ�������ָ��Ϳ�����
	//�����ǿ��ԸĽ��ĵط�
	Filter dFilter;
	SFilter dFilter1;
	//�ѿ�������
	int type;	
	//�������ǲ����Ѿ������
	bool IsDeal;
	//�Ƿ����SVMģ��
	bool IsOPSVMModel;
	//�Ƿ����У��ģ�ͽ�ģ����
	bool IsOpCData;
	//�ǲ������MS1���ݣ������ʽ��ʲô
	int ms1format;
	//MS2�����ʽ�����������
	int ms2format;
	//�Ƿ�ʹ��FT״̬����
	bool IsUseStatsPar;
	//У��ģ��ѡ�񣬰��������κ�ģ�ͣ�����ģ�ͺ�svmģ��3��
	int ModelType;
public:
	ThreadParm();
	ThreadParm(const ThreadParm &r);
	~ThreadParm();	
	//����/�ͷű�������Ĵ洢�ռ�
	void AllocRepBuff(int n);
	void FreeRepBuff();
	//��д/��ȡ��������
	bool SetRepBuff(int idx,CReport *rt);
	bool GetRepBuff(int idx,CReport *rt);
	//���ñ���״̬
	void ResetRepState(int idx);
	ThreadParm &operator=(const ThreadParm &r);
};

//�߼���������
class AdvPars
{
public:
	int Min_ISO;//��Сͬλ�ط���Ŀ
	double Min_Rel_Int;//��С����ź�ǿ��
	double Iso_GD;//��С����Ŷ�
	int Xic_Int_Max;//XIC���������gap�����������ô����ͽض�
	double RT_Max_Int;//XIC�����������ʱ����
	double Ch_RT_Min;//���ȷ����������Сֵ
	double Ch_RT_Max;//���ȷ����������Сֵ
	double Pre_dm;//Ԥ����ĸ������������XIC������ʹ��
	//double IsoDiff[MAX_ISO];
	//double stats_Iso_MET[MAX_ISO];
	double SVM_C;//svmѵ��������Խ����Ҫʱ��Խ��
	//double SVM_gama;
	bool PIonResearch;//�Ƿ����ĸ������չ����
	bool XICCalibrate;//�Ƿ�ִ��XICУ��
	long TrainSize;//svmѵ�����������������������������ڸ���ֵ���������ѡ��
	int Max_Ch;//���Ŀ��ǵ��
	int MinIsoNum;//��չ����ĸ����ʱ��ͬλ�ط���Сֵ
	double RIntCTH;//��չ����ʱ������ź�ǿ����Сֵ
	double GdCtH;
	int LExt;
	int RExt;
	//to limit the search work
	double MinRT;
	double MaxRT;
	//
	double Min_SN;//��������
	//��̬XIC�ضϲ���
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
