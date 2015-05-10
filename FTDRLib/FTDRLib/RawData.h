#pragma once
#include "ThermoInterface.h"
#include "SVMModel.h"
#include "CalScan.h"
#include "commonddef.h"
#include "HistList.h"
#include "PepXML.h"
#include "LinModel.h"
#include "dSGSmooth.h"

#define STD_TUNE 2.698
#define STD_EPS 1e-16;
#define ERROR_MASS -1.0

//#define CASE_TEST 1
//#define DEBUG_CASE 1

class MS2Header
{
public:
	double pmz;
	long scan;
	double IW;
	int ch;
	long MS1Scan;
public:
	MS2Header();
	~MS2Header();
	MS2Header(const MS2Header &cp);
	MS2Header &operator =(const MS2Header &cp);
	//must be call by specfic function, not generally
	void Write2File(FILE *fp);
	void ReadFromFile(FILE *fp);
};

class CalReturn
{
public:
	double cal_mean;
	double cal_std;
	int charge;
	double rel_int;
	double iso_goodness;
	int iso_num;
	int spectrum_num;
	double abs_int;
public:
	CalReturn();
	~CalReturn();
	CalReturn(const CalReturn &cp);
	CalReturn &operator=(const CalReturn &cp);
	CalReturn &operator=(const CalData &cp);

};

class RawData
{
public:
	vector<size_t> ColParIdx;
	size_t TailerIdx[2];
	vector<size_t> SelParIdx;
	vector<string> ParNames;
private:
	int Error_code;
	int CalTotal;
	int modeTotal;
private:
	Filter cFilter;	
	SFilter cFilter1;
//for calibrate, the raw file can be opened for one time and used in the fowllowing calculation
//if TIF==NULL, the interface is not initial, else we should delete the TIF object in the de-construct function
public:
	myHistogram mtA;//the mass error distribution before calibration
	myHistogram mtB;//the mass error distribution after calibration
	double METModel[4];
	bool IsMETModel;
	void GetPreHis(vector<CalData> &XICList);
	void GetCalHis(vector<CalData> &XICList);
private:
	ThermoInterface *TIF;
	void InitialCalOutput(char *raw);
public://for output debug data
	FILE *fp_cal_data;
	bool IsOutput;
private:
	vector<CalScan> MSList;
	vector<MS2Header> MS2Scans;
	//double current_mz;
	//int current_ch;	
	char raw_name[MAX_PATH];
	double IsoDisT[MAX_ISO];
	SVMModel *CalModel;
	MSInstrumentModelType instrument_T;
private:
	//find the adjucent ms 1 scan from ms2 scan
	long seekScan(long scannum,size_t a, size_t b);
	vector<long> ScanMap;	
public:
	char IntRModelData[MAX_PATH];	
	bool OIntR_Data;
public:
	LinModel LModel;
	int ModelType;
public:
	RawData(void);
	~RawData(void);		
	int ProcessOneXICBack(CalData &ct,int MS1Scan,double bReturn[2]);
	int ProcessOneXIC_Mono(CalData &ct,int MS1Scan,double bReturn[2]);
	int ProcessOneXIC(CalData &ct,int MS1Scan,double bReturn[2]);
	int ProcessOneXIC_Lin(CalData &ct,int MS1Scan,double bReturn[2]);	
	int ProcessOneXIC_Non(CalData &ct,int MS1Scan,double bReturn[2]);	

	bool XICCalibrate(size_t MS2,vector<CalReturn>& bReturn);	
	bool XICCalibrate_Lin(size_t MS2,vector<CalReturn>& bReturn);	
	bool XICCalibrate_Mono(size_t MS2,vector<CalReturn>& bReturn);

	bool XICCalibrate_Non(size_t MS2,vector<CalReturn>& bReturn);
	bool XICCalibrateC_Non(size_t MS2,vector<CalReturn>& bReturn);
	bool XICCalibrateS_Non(size_t MS2,vector<CalReturn>& bReturn);

	bool XICCalibrateC(size_t MS2,vector<CalReturn>& bReturn);
	bool XICCalibrateS(size_t MS2,vector<CalReturn>& bReturn);

	bool XICCalibrateC_Lin(size_t MS2,vector<CalReturn>& bReturn);
	bool XICCalibrateS_Lin(size_t MS2,vector<CalReturn>& bReturn);

	bool SCalibrate(size_t MS2,vector<CalReturn>& bReturn);
	bool SCalibrateC(size_t MS2,vector<CalReturn>& bReturn);
	bool SCalibrateS(size_t MS2,vector<CalReturn>& bReturn);	

	bool SCalibrate_Lin(size_t MS2,vector<CalReturn>& bReturn);
	bool SCalibrateC_Lin(size_t MS2,vector<CalReturn>& bReturn);
	bool SCalibrateS_Lin(size_t MS2,vector<CalReturn>& bReturn);

	double ECalibrate(size_t MS2);
	double ECalibrate_Lin(size_t MS2);
	void GetXICData(vector<CalData> &XICList,double pmz,int ch,int MS2Scan);
	void GetData(vector<CalData> &SList,double pmz,int ch,int MS2Scan);
	//adding version
	void GetXICData(vector<CalData> &XICList,double mze,int ch,int MS2Scan,double mzt);
	void GetXICData_dMET(vector<CalData> &XICList,double mze,int ch,int MS2Scan,double mzt);
	void GetXICData_dMET_NoLim(vector<CalData> &XICList,double mze,int ch,int MS2Scan,double mzt);
	void GetData_dMET(vector<CalData> &SList,double mze,int ch,int MS2Scan,double mzt);
	bool PreCalData(CalData &ct,double mze,int ch,int MS2Scan,double mzt);
	double GetMaxInt(vector<CalData> &XICList);
	void GetPossChFromPmz(double pmz,int scan,vector<int> &CHL);
	void SortMerge(vector<mPeak> &tmpPKL);
	double GetSampleDm(double mass);
	bool ModelBuilding(vector<CalData> &XICList);	
	bool ModelBuilding_Non(vector<Pre_calData> &IDCalData);
	bool ModelBuilding_dMET(vector<Pre_calData> &IDCalData);
	bool ModelBuilding_Lin(vector<Pre_calData> &IDCalData);
	bool ModelBuilding_X(vector<Pre_calData> &IDCalData);
	bool ModelBuilding_S(vector<Pre_calData> &IDCalData);
	bool ModelBuilding_LinX(vector<Pre_calData> &IDCalData);
	bool ModelBuilding_LinS(vector<Pre_calData> &IDCalData);	
	//for limit the SVM training data, too many will be error
	void RandomSelect(vector<CalData> &XICList);
	bool EstimateCalMET(vector<Pre_calData> &IDCalData);
	bool EstimateCalMET_Mono(vector<Pre_calData> &IDCalData);
	bool EstimateCalMET_dMET(vector<Pre_calData> &IDCalData);	
	bool EstimateCalMET_Lin(vector<Pre_calData> &IDCalData);
	bool EstimateCalMET_S(vector<Pre_calData> &IDCalData);
	bool EstimateCalMET_LinS(vector<Pre_calData> &IDCalData);
	void GetCalHis(vector<CalData> &XICList,double bT[2]);
	double GetCalPmz(vector<CalData> &XICList);
	bool GetCalPmz(vector<CalData> &XICList,double bRT[2]);
	double GetCalPmz_Lin(vector<CalData> &XICList);
	double GetCalPmz_Non(vector<CalData> &XICList);
	bool GetCalPmz_Lin(vector<CalData> &XICList,double bRT[2]);
	bool GetCalPmz_Non(vector<CalData> &XICList,double bRT[2]);
public:
	double MET;	
	double ppmLV,ppmHV;
	void GetPreMET(vector<Pre_calData> &IDCalData);
	bool Calibrate(char *mgfname);
	bool Calibrate_Lin(char *mgfname);
	bool Calibrate_Non(char *mgfname);
	bool Calibrate_SVM(char *mgfname);
	bool Calibrate_Complex(char *mgfname);
	bool Calibrate_Simple(char *mgfname);
	bool Calibrate_Com_Lin(char *mgfname);
	bool Calibrate_Sim_Lin(char *mgfname);
	bool Calibrate_Mono(char *mgfname);
	bool Calibrate(char *rawfile,vector<Pre_calData> &IDCalData,char *mgfname);	
public://debug functions
	void EnableOutput(char *fname);
	void CloseOutput();
	int Loaddata(char *sraw);
	long GetScanNum(CString sFName);
	int CollectDataS(CString RPath,vector<Pre_calData> &IDCalData);
	int CollectDataP(CString RFile,vector<Pre_calData> &IDCalData);
	int CollectDataM(CString RFile,vector<Pre_calData> &IDCalData);
	int CollectDataX(CString RFile,vector<Pre_calData> &IDCalData);
	bool Initial(char *rawfname);
public:
	int GetErrorCode();
	void SetErrorCode(int errors);
	void GenReport(CReport *rep);
	void InitialFilter(Filter *dFilter,SFilter *dFilter1);
	void OutputModel(char *svmmodel);
public:
	//for debug to op the SVM parameters
	bool OutPutCalData(vector<Pre_calData> &IDCalData,char *fname);
	void OutPutCalData(vector<CalData> &XICList,FILE *fp);
	void GeneratePurF(char tmpPureRawName[],char raw_names[]);
public:
	//to estimate a intensity relative model for MET
	void SetOutputRData(char *fname);
public:
	//to clear for a new states
	void clear();
	SVMModel *GetSVMPointer();
public:	
	void PreSuvey();
	void featureSel(vector<CalData> &XICList);
	bool IsUseExtendPar;	
public:
	int ppbLevel;
public:
	vector<string> DInfBuf;
	void OutPutInf();
public:
	double MzIntModel[4];
	double GetSmDm(double mass);
#ifdef CASE_TEST
	FILE *ParentIons;
#endif

public:
	bool SaveRawData(char *fname);
	bool LoadRawData(char *fname);
public:
	int ms2format;
	bool OutputMS2(FILE *MGFfp,FILE *IonsFile,size_t MS2,vector<CalReturn> &bReturn,char PurerawName[]);
	bool IsParentExist(vector<CalReturn> &bReturn,double pmz);
public:
	int ms1format;
	void GetPurRawName(string &purname);
	bool OutPutMS1(string outpath);
public:
	DiffCut DFT;
	void InitialDFTModel();
};
