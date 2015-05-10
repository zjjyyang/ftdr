// mPeak.cpp: implementation of the mPeak class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mPeak.h"
#include "math.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mPeak::mPeak()
{
	dMass=0;
	dIntensity=0;
	rInt=0;
	baseLine=0;
	noise=0;	
}

mPeak::~mPeak()
{

}

mPeak &mPeak::operator=(const mPeak &r)
{
	dMass=r.dMass;
	dIntensity=r.dIntensity;
	rInt=r.rInt;
	baseLine=r.baseLine;
	noise=r.noise;
	return *this;
}

mPeak::mPeak(const mPeak &r)
{
	dMass=r.dMass;
	dIntensity=r.dIntensity;
	rInt=r.rInt;
	baseLine=r.baseLine;
	noise=r.noise;	
}

void mPeak::Write2File(FILE *fp)
{
	fwrite(&dMass,sizeof(double),1,fp);
	fwrite(&dIntensity,sizeof(double),1,fp);
	fwrite(&rInt,sizeof(double),1,fp);
	fwrite(&baseLine,sizeof(double),1,fp);
	fwrite(&noise,sizeof(double),1,fp);
}

void mPeak::ReadFromFile(FILE *fp)
{
	fread(&dMass,sizeof(double),1,fp);
	fread(&dIntensity,sizeof(double),1,fp);
	fread(&rInt,sizeof(double),1,fp);
	fread(&baseLine,sizeof(double),1,fp);
	fread(&noise,sizeof(double),1,fp);
}