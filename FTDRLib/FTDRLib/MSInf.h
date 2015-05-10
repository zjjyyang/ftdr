#pragma once
#include "commonfit.h"
class MSInf
{	
public:
	bool isCentroided;
	DataPeak *pkl;
	int PNum;
	long ScanNumber;
	double RT;
	MSInf(void);
	~MSInf(void);
	MSInf(const MSInf &r);
	MSInf & operator =(const MSInf &r);
	bool SetPKL(double *pdval,int dim);
	bool SetPKL(DataPeak *tmpPKL,int count);
	bool SetPKL(vector<double> &mz,vector<double> &abu);
	double GetBasePeak();
	int Centriod(string instrument);
	bool OutPutData(FILE *fp);
	double Signal_Th;
	double RSlower;
	double GetSampleDm(string instrument,double mass);
	int FindByPmz(FILE *fp);
	////////////////////
	int Min_Iso_Num;
	double GD_CUT;
	double MassErrorTol;
	void RemUsed(vector<Peak> &tmpPKL);
	int FindIsoCPmz(double BaseInt,FILE *fp,vector<Peak> &tmpPKL);
};
