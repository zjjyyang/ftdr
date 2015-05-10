// iPair.cpp: implementation of the iPair class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "iPair.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

iPair::iPair()
{
	RT=0;
	dm=0;
	ScanNum=0l;
}

iPair::~iPair()
{

}

iPair &iPair::operator=(iPair &r)
{
	for(int i=0;i<3;i++)
	{
		LPeak[i]=r.LPeak[i];
		HPeak[i]=r.HPeak[i];
	}
	RT=r.RT;
	dm=r.dm;
	ScanNum=r.ScanNum;
	return *this;
}

double iPair::GetRatioHL(double r)
{
	double rh=HPeak[0].dIntensity-HPeak[0].dBaseLine-HPeak[0].dNoise;
	double rl=LPeak[0].dIntensity-LPeak[0].dBaseLine-LPeak[0].dNoise;
	if(rl<=0) return 0;
	rh=rh-(LPeak[2].dIntensity-LPeak[2].dBaseLine-LPeak[2].dNoise)*r;
	if(rh<0) return 0;
	return rh/rl;
}

double iPair::GetSN()
{
	if(LPeak[0].dNoise<=0) return 0;
	return LPeak[0].dIntensity/LPeak[0].dNoise;
}