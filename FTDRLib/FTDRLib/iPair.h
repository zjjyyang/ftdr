// iPair.h: interface for the iPair class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_IPAIR_H__7F4F7311_A748_4C9D_9CD7_99FF77946C7A__INCLUDED_)
#define AFX_IPAIR_H__7F4F7311_A748_4C9D_9CD7_99FF77946C7A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "mPeak.h"

class iPair  
{
public:
	mPeak LPeak[3];
	mPeak HPeak[3];	
	double RT;
	double dm;
	long ScanNum;	
	iPair();
	~iPair();
	iPair &operator=(iPair &r);
	double GetRatioHL(double r);
	double GetSN();
};

#endif // !defined(AFX_IPAIR_H__7F4F7311_A748_4C9D_9CD7_99FF77946C7A__INCLUDED_)
