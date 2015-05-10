// iProtein.h: interface for the iProtein class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_IPROTEIN_H__1948451F_E4C6_4F89_9AAA_453C75A8F4C2__INCLUDED_)
#define AFX_IPROTEIN_H__1948451F_E4C6_4F89_9AAA_453C75A8F4C2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "iPeptide.h"
#include "afxtempl.h"

class iProtein  
{
private:
	iPeptide *ipl;
	int count;
public:  
	double ratio;
	double std;
	void UpdateScore();
	float tscore;
	float mass;
	int pcount;
	CString pAcc;
	iProtein();
	iProtein(int num);
	iProtein(int num,iPeptide *pl);
	iProtein(CList<iPeptide,iPeptide&> *pl);
	iProtein(CList<iPeptide,iPeptide&> &pl);
	int GetCount();
	bool SetPL(int num,iPeptide *pl);
	bool SetPL(CList<iPeptide,iPeptide&> *pl);
	bool SetPL(CList<iPeptide,iPeptide&> &pl);
	bool GetPL(CList<iPeptide,iPeptide&> *pl);
	bool GetPL(CList<iPeptide,iPeptide&> &pl);
	int GetPL(int num,iPeptide *pl);
    bool SetAt(int idx,iPeptide &r);
	iPeptide GetAt(int idx);
	bool GetAt(int idx,iPeptide &r);
	iPeptide operator [](int idx);
	iProtein &operator =(iProtein &r);
	~iProtein();
};

#endif // !defined(AFX_IPROTEIN_H__1948451F_E4C6_4F89_9AAA_453C75A8F4C2__INCLUDED_)
