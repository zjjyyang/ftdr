// OutRecord.h: interface for the OutRecord class.

#pragma once

#include "SRecord.h"

#define MAX_LINEW 1000

class SFilter
{
public:
	SFilter();
	~SFilter();
	double XCorr[3];
	double detCn;
	int RSp;
	//int MinPepLen;
	//int MaxPepLen;
	SFilter(const SFilter &r);
	SFilter &operator=(const SFilter &r);
	bool IsPassFilter(SRecord *pep,int ch);
	void WriteToFile(CStdioFile &dFile);
	bool LoadFromFile(CStdioFile &dFile);
	bool IsValidate();
};

class OutRecord  
{
private:
	SRecord *SR;
	int count;
public:
	float GetdetCn();
	CString sFName;
	int MatchPep;
	float TotalInt;
	float EMH;
	float LowSp;
	int Charge;
	OutRecord();
	~OutRecord();
	int GetCount();
	bool GetAt(int idx,SRecord *rec);
	bool GetAt(int idx,SRecord &rec);
	bool SetAt(int idx,SRecord *rec);
	OutRecord& operator=(OutRecord &r);
	bool ReadFromFile(CString fname);
	bool RemTailSP(char *str);
	bool StrTailFind(char *str,char *substr);
	bool ReadFromFile(FILE *fp,char *outtitle);
};

