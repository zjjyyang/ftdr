#pragma once
#include "mPeak.h"
#include "commonddef.h"
#include "Scan.h"
#include "IsoGroupD.h"
class CalScan
{
public:	
	vector<double>StatusPar;
	double ISt;
	double ESt;
	long ScanNum;
	double RT;	
	double BaseInt;
	int GetParIdx(string &sLabel);
	int TFindH(int b,int e,double fmz);
	int TFindL(int b,int e,double fmz);	
	int LocatePmz(double DM,double pmz);
	int LocatePmz(double begin,double end,double pmz);		
	bool Findmz(mPeak &mt,double pmz,double dm);	
	bool FindISO(CalData &it,int ch,double pmz);
	bool FindISO_dMET(CalData &it,int ch,double pmz,double ppmL,double ppmH);
	//support to find the overlap peaks
	bool ExtendFindISO(vector<CalData> &PIso,double pmz,double IW);		
	bool FindALLISO(vector<CalData> &PIso);
	void RemUsed(vector<sPeak> &tmpPKL);
	double GetTIC();
	void PackageIt(CalData &it);
	int FindIsoCPmz(vector<sPeak> &tmpPKL,vector<CalData> &PIso);
public:	
	double Iso_win;
	double Pre_dm;
	int Min_Iso_Num;
	vector<mPeak> PKL;
	CalScan(void);
	~CalScan(void);
	CalScan(const CalScan &cp);
	CalScan &operator=(const CalScan &cp);
	bool Convert(Scan *cs,vector<size_t> &PreSuvList);	
	bool Convert(Scan *cs,vector<size_t> &StatusIdx,size_t TailerIDX[2]);	
	void RetrivalPeaks(double min_mz,double max_mz,vector<mPeak> &tmpPKL);
	void InitialSTSisze(size_t sg);
public:
	fpos_t Write2File(FILE *fp);
	void ReadFromFile(FILE *fp);
};

