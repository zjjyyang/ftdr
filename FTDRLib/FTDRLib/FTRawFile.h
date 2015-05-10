// FTRawFile.h: interface for the FTRawFile class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

//#include "afxtempl.h"
#include "mPeak.h"
#include "xrawfile2.h"
//#include "commonddef.h"
//#include <vector>
//using std::vector;
#define FT 1
#define ORBITRAP 2

#define STD_TUNE 2.698
#define STD_EPS 1e-16;

#define ERROR_MASS -1.0

class FTRawFile : public IXRawfile  
{
private://for debug
	FILE *fp_cal_data;
	bool IsOutput;
private:
	long SNLast,SNFirst;
	double bMass;
	int charge;	
	vector<mPeak> MSSpec;
	char sFName[MAX_PATH];
	double IsoDisT[MAX_ISO];

	double GetMSTIC();
	bool RefinePISO(IsoCluster &it);
	bool PreRefinePISO(IsoCluster &it);
	double PreFindMIntSim();
	void CalDefIsoDis(double mass,double IsoAreaT[6]);
	bool FindMInt(IsoCluster &it);	
	double PreFindMInt(IsoCluster &it);
	bool XICCalibrate(double C[5],FILE *MGFfp,long ScanNum);
	bool XICCalibrateExt(double C[pdNum],FILE *MGFfp,long ScanNum);
	double XICCalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI);
	bool XICCalibrateBack(double C[5],FILE *MGFfp,long ScanNum);

	bool SCalibrate(double C[5],FILE *MGFfp,long ScanNum);
	bool SCalibrateExt(double C[pdNum],FILE *MGFfp,long ScanNum);
	double SCalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI);
	bool ECalibrate(double C[5],FILE *MGFfp,long ScanNum);
	bool ECalibrateExt(double C[pdNum],FILE *MGFfp,long ScanNum);
	double ECalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI);
	void GetXICDataExt(vector<IsoCluster> *XICList,vector<double> *dMass,double C[pdNum]);
	void GetXICData(vector<IsoCluster> *XICList,vector<double> *dMass,double C[5]);	
	void GetIsoCData(IsoCluster *it,vector<double> *dMass,double C[5]);
	void GetIsoCDataExt(IsoCluster *it,vector<double> *dMass,double C[pdNum]);	
	double match(mPeak ISOP[MAX_ISO],double IsoDisT[MAX_ISO]);
	void CalMST(vector<double> *dMass,double *m,double *s);
	void rCalMST(vector<double> *dMass,double *m,double *s);
	int TFindL(vector<mPeak> &pkl,int b,int e,double fmz);
	bool ReadMS(long ScanNum,vector<mPeak> &MSSpecTmp);
	bool ReadMSLabel(long ScanNum,vector<mPeak> &MSSpecTmp);
	double GetPMassL(long ScanNum,double mzeP);
	int MSLevel(long ScanNum);
	bool GetMS2pMass(long ScanNum);		
	//added on 2010.6.14 by zhangjy to find the best parent mass 
	double PreFindMIntWt();
	double PreFindMIntFt();
	double PreFindMIntMP();
	bool FindMIntFit(IsoCluster &it);

	bool PreFindMIntFtISO(IsoCluster &it);
	bool FindMIntFtISO(IsoCluster &it);

	bool FindMass(double massE[5],long ScanNum,double massP,int ch);
	double GetSampleDm(int type,double mass);
	int FindPKPos(int begin, int end,vector<mPeak> PPos);

	int LocatePmass(double begin,double end,double pmass);
	int LocatePmassH(double begin,double end,double pmass);
	int LocatePmass(double begin,double end,double pmass,double DM);

	int TFindL(int b,int e,double fmz);
	int TFindH(int b,int e,double fmz);

	int PreFindMIntFt(mPeak &mt);
	int FindMIntFt(mPeak &mt);
	void GetCOS(IsoCluster &it,int num);

	bool PreFindMIntFtNew();

	double w_mean(vector<double> *x, char *W);
	double w_std(vector<double> *x,char *W,double m);
	void UpdateResW(vector<double> *x,char *W,double mean,double sigma);
	bool GetFTStatusData(long ScanNum,CalData &dt);
	int GetParIdx(CString &sLabel);
	int GetParIdxAll(CString &sLabel);
	bool GetFTEJTData(long ScanNum,CalData &dt);

	bool GetFTStatusData(long ScanNum,IsoCluster &dt);	
	bool GetFTEJTData(long ScanNum,IsoCluster &dt);
	double BasePeakInt();
	//added on 2010.11.10
	void CalSKPar(vector<double> x,vector<double>y,double beta[4],double ReturnPar[4]);
	void CalSKPar(vector<double> x,vector<double>y,vector<double> beta,double ReturnPar[4]);
public:
	double MET;
	FTRawFile();
	~FTRawFile();
	bool Initial(CString sRawfile);
	bool GetTIC(long ScanNum,double *TIC, double *RT);
	int Calibrate(double C[5],char *MGFfile);
	int CalibrateExt(double C[pdNum],char *MGFfile);	
	double CalibrateExt(double C[pdNum],long ScanNum,double mzinitial,int ch);
	double GetGCError(double C[5], long scan,double pmass,int ch);
	double GetGCErrorExt(double C[pdNum], long ScanNum,double pMZ,int ch);
	int SimpleCalibrate(double C[5],char *MGFfile);	
	int SimpleCalibrateExt(double C[pdNum],char *MGFfile);
	double SimpleCalibrateExt(double C[pdNum],long ScanNum,int ch,double mzI);
	void EnableOutput(char *fname);
	void CloseOutput();
	bool GetFTingData(double pMass,long ScanNum,CalData &dt);
	bool GetFTingData(long ScanNum,CalData &dt);
	bool GetFTingData(double pMass,long ScanNum,CalData &dt,int ch);
	bool GetFTExtraData(long ScanNum,CalData &dt);//for the MS spectrum	

	//for raw calibration
	bool GetFTingData(double pMass,long ScanNum,IsoCluster &dt);
	bool GetFTingData(long ScanNum,IsoCluster &dt);//use bMass
	bool GetFTingData(double pMass,long ScanNum,IsoCluster &dt,int ch);
	bool GetFTExtraData(long ScanNum,IsoCluster &dt);//for the MS spectrum	
	
	//to investigate the peak skew
	//added on 2010.11.10
	bool GetSkewData(long ScanNum,double pMZ);

	//added to load all parameters, on 2010.12.3
	bool GetFTingDataAll(double pMass,long ScanNum,CalData &dt,int ch);
	bool GetFTStatusDataAll(long ScanNum,CalData &dt);
	bool GetFTingDataAll(double pMass,long ScanNum,CalData &dt);
};
