#if !defined(AFX_MSQADOC_H__6AB2CB25_5689_4E4D_957D_EA88FA4493D3__INCLUDED_)
#define AFX_MSQADOC_H__6AB2CB25_5689_4E4D_957D_EA88FA4493D3__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "commonfit.h"
#include <vector>
using std::vector;
#include "afxtempl.h"
#include "MSInf.h"


class IdxRaw
{
private:
	vector<XICs> _pXICList;
	vector<Peak> tmpPKL;
	size_t CurrentMS;
	double CurrentBestGD;
	void GetSNBL(long nScanNumber);	
	int FindIsoCNew(double BaseInt);
	void FindMInt(DataPeak *pkl,int pnum);
	double GetBasePeak(DataPeak *pkl,int pnum);	
	int FilterPP(XICs &xct,vector<XICs> *XICList);
	void AddXICPeaks(XICPeak &xt,int CH);
	void AddXICPeaks(XICPeak &xt);
	void RemUsed();
	int FindByPmz(double pmz,DataPeak *pkl,int pnum,FILE *fp);
	int TFindL(DataPeak *pkl,int b,int e,double fmz);
	int IDXFromScanNum(long ScanNum);
	int TFindL(int b,int e,long scan);
	int FindIsoCPmz(double BaseInt,FILE *fp);
	int FindIsoCPmz();
	int FindIsoCPmz(int *PreData,int PreNum,FILE *fpout);
public:	
	IdxRaw(void);
	~IdxRaw(void);	
	vector<MSInf> MSList;
	double MinS_N;
	double MassErrorTol;
	size_t Min_Iso_Num;
	size_t Min_Peaks_Num;
	double GD_CUT;
	double ISO_Win;
	int FindPP(vector<XICs> *XICList);	
	int ISObyPmz(long ScanNum,double pmz,int ch,FILE *fp);
	bool ISObySingleMS(char *MSFile);
	bool ISObySingleMSNew(char *MSFile,char *outfile);
	void Close();	
};
#endif // !defined(AFX_SMVIEW_H__5D0F1C35_8F84_4263_AD91_8F2CE9E14388__INCLUDED_)