// iPeptide.cpp: implementation of the iPeptide class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "iPeptide.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


iPeptide::iPeptide()
{
	 Observed=0;
	 CalcMass=0;
	 ObsrvMass=0;
	 dm=0;
	 miss=0;
	 score=0;
	 rank=0;
	 color=COLOR_NON;
	 IsBold=false;
	 IsChecked=false;	
	 Expect=1;
	 ScanNumber=0;
}

iPeptide::~iPeptide()
{

}

iPeptide &iPeptide::operator=(iPeptide &r)
{
	 Observed=r.Observed;
	 CalcMass=r.CalcMass;
	 ObsrvMass=r.ObsrvMass;
	 Expect=r.Expect;
	 dm=r.dm;
	 miss=r.miss;
	 score=r.score;
	 rank=r.rank;
	 color=r.color;
	 IsBold=r.IsBold;
	 IsChecked=r.IsChecked;	
	 ScanNumber=r.ScanNumber;
	 fName=r.fName;
	 Sequence=r.Sequence;
	 Modif=r.Modif;
	 return *this;
}

void iPeptide::Clear()
{
	 Observed=0;
	 CalcMass=0;
	 ObsrvMass=0;
	 Expect=0;
	 dm=0;
	 miss=0;
	 score=0;
	 rank=0;
	 color=COLOR_NON;
	 IsBold=false;
	 //IsChecked=false;
	 //ScanNumber=0;
	 //fName.Empty();
	 Sequence.Empty();
	 Modif.Empty();
}
