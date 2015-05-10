// alQuant.h: interface for the CalQuant class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include "RawData.h"

#define WM_DISPLAY (WM_USER+101)

//the main task is to load the identifications
class CalQuant  
{
private:	
	RawData *RWCal;	
	void GenMF(ThreadParm *parm,CString &CalModelDataFile);
	void GetFNameC(CString fpath, CString &fname);
	void GetFName(CString fpath,CString &fname);
	void GetFNameF(CString fpath, CString &fname);
	void GetFNameMd(CString fpath, CString &fname);	
public:	
	UINT ThreadID;		
	CalQuant(UINT ID);
	~CalQuant();
	bool IsOutSuccess;
	void ExecuteCal(ThreadParm *parm);
	void ExecuteCal_Lin(ThreadParm *parm);
	void ExecuteCal_SVM(ThreadParm *parm);
	void ExecuteCal_Non(ThreadParm *parm);
	void ExecuteCal_Mono(ThreadParm *parm);
	void PostFalseMsg(WPARAM w_Report);
	void GenReport(CReport *rep);
	void SaveSVMModel(char *fname);
	void OutPutCalData(ThreadParm *parm);
};

