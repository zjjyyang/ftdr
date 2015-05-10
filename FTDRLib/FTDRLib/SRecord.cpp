// SRecord.cpp: implementation of the SRecord class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRecord.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SRecord::SRecord()
{
	RankSp=0;
	MH=0;
	detCn=0;
	XCorr=0;
	Sp=0;
	mIons=0;
	tIons=0;
}

SRecord::~SRecord()
{

}

SRecord &SRecord::operator =(SRecord &r)
{
	RankSp=r.RankSp;
	MH=r.MH;
	detCn=r.detCn;
	XCorr=r.XCorr;
	Sp=r.Sp;
	mIons=r.mIons;
	tIons=r.tIons;
	Ref=r.Ref;
	Seq=r.Seq;
	return *this;
}
