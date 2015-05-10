// MResult.h: interface for the MResult class.

#pragma once

#include "iPeptide.h"
#include "afxtempl.h"

#define LIGHT_RED 0x01
#define BOLD_BLACK 0x02
#define LIGHT_BLACK 0x04

class Filter
{
public:
	Filter();
	~Filter();
	double DefaultMin;
	double MinPepScore;
	int MaxRank;
	int MaxMiss;
	double MET;
	BOOL IsUseDefault;
	Filter(const Filter &r);
	Filter &operator=(const Filter &r);
	bool IsPassFilter(iPeptide *pep);
	void WriteToFile(CStdioFile &dFile);
	bool LoadFromFile(CStdioFile &dFile);
	bool IsValidate();
};

class MResult  
{
private:
	iPeptide *PEPLIST;
	int count;
	bool ReadOnePep(CString sTemp,iPeptide &ptemp);
	bool Isppm;
public:
	bool LoadData(CString fname);
	int GetVersion(CString fname);
	Filter cFilter;
	bool IsCarb;
	void RemoveMarker(CString &sTemp);
	int GetScanNumber(CString item);
	void Clear();
	bool RemoveRon(CList<iPeptide,iPeptide&> *ipl);
	int GetCount();
	MResult();
	MResult(int num);
	bool SetPL(int num,iPeptide *pl);
	bool SetPL(CList<iPeptide,iPeptide&> *pl);
	bool SetPL(CList<iPeptide,iPeptide&> &pl);
	bool GetPL(CList<iPeptide,iPeptide&> *pl);
	bool GetPL(CList<iPeptide,iPeptide&> &pl);
	int GetPL(int num,iPeptide *pl);
    bool SetAt(int idx,iPeptide &r);
	bool GetAt(int idx,iPeptide &r);
	bool SetAt(int idx,iPeptide *r);
	bool GetAt(int idx,iPeptide *r);
	MResult &operator =(MResult &r);
	int ReadFromFile(LPCTSTR fname);
	int ReadFromFileNew(LPCTSTR fname);
	bool ReadOnePep2(CString sTemp,iPeptide &ptemp);
	int ReadFromFile2(LPCTSTR fname);
	void RemoveSPs(CString &sTemp);
	void ReplaceTokens(CString &sTemp);
	~MResult();

	//new common methods to read the html file of mascot
	bool ReadAtable(CStdioFile &dFile,CString &str);
	int SubStrCount(CString &str,char *substr);
	void MakeLower(CString &str);
	int LoadFromDat(LPCTSTR fname);
	bool ValidateEXT(LPCTSTR fname,const char *EXT);
};

