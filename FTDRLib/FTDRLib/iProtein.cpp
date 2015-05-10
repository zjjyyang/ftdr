// iProtein.cpp: implementation of the iProtein class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "iProtein.h"
#include "math.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

iProtein::iProtein()
{
	ipl=NULL;
	ratio=0;
	std=0;
	count=0;  
}

iProtein::~iProtein()
{
	if(count>0&&ipl!=NULL) delete []ipl;
}


iProtein::iProtein(int num)
{
	if(num>0) 
	{
		ipl=new iPeptide[num];
		count=num;
	}
	else
	{
		ipl=NULL;
		count=0;  
	}
}

iProtein::iProtein(int num,iPeptide *pl)
{
	if(num>0) 
	{
		ipl=new iPeptide[num];
		count=num;
		for(int i=0;i<num;i++) ipl[i]=pl[i];
		UpdateScore();
	}
	else
	{
		ipl=NULL;
		count=0;  
	}
}

iProtein::iProtein(CList<iPeptide,iPeptide&> *pl)
{
	int num=pl->GetCount();
	if(num>0) 
	{
		ipl=new iPeptide[num];
		count=num;
		POSITION pos=pl->GetHeadPosition();
		for(int i=0;i<num;i++) ipl[i]=pl->GetNext(pos);
		UpdateScore();
	}
	else
	{
		ipl=NULL;
		count=0;  
	}
}

iProtein::iProtein(CList<iPeptide,iPeptide&> &pl)
{
	int num=pl.GetCount();
	if(num>0) 
	{
		ipl=new iPeptide[num];
		count=num;
		POSITION pos=pl.GetHeadPosition();
		for(int i=0;i<num;i++) ipl[i]=pl.GetNext(pos);
		UpdateScore();
	}
	else
	{
		ipl=NULL;
		count=0;  
	}
}

bool iProtein::SetPL(int num,iPeptide *pl)
{
	if(count>0&&ipl!=NULL)
	{
		delete []ipl;
		count=0;
		ipl=NULL;
	}
	if(pl==NULL) return false;
	if(num>0) 
	{
		ipl=new iPeptide[num];
		count=num;
		for(int i=0;i<num;i++) ipl[i]=pl[i];
		UpdateScore();
		return true;
	}
	else return false;
}

bool iProtein::SetPL(CList<iPeptide,iPeptide&> *pl)
{
	if(count>0&&ipl!=NULL)
	{
		delete []ipl;
		count=0;
		ipl=NULL;
	}
	if(pl==NULL) return false;
	int num=pl->GetCount();
	if(num>0) 
	{
		ipl=new iPeptide[num];
		count=num;
		POSITION pos=pl->GetHeadPosition();
		for(int i=0;i<num;i++) ipl[i]=pl->GetNext(pos);
		UpdateScore();
		return true;
	}
	else return false;	
	
}

bool iProtein::SetPL(CList<iPeptide,iPeptide&> &pl)
{
	if(count>0&&ipl!=NULL)
	{
		delete []ipl;
		count=0;
		ipl=NULL;
	}
	int num=pl.GetCount();
	if(num>0) 
	{
		ipl=new iPeptide[num];
		count=num;
		POSITION pos=pl.GetHeadPosition();
		for(int i=0;i<num;i++) ipl[i]=pl.GetNext(pos);
		UpdateScore();
		return true;
	}
	else return false;

}

bool iProtein::GetPL(CList<iPeptide,iPeptide&> *pl)
{
	if(pl==NULL) return false;
	pl->RemoveAll();
	for(int i=0;i<count;i++)
	{
		pl->AddTail(ipl[i]);
	}
	return true;
}

bool iProtein::GetPL(CList<iPeptide,iPeptide&> &pl)
{
	pl.RemoveAll();
	for(int i=0;i<count;i++)
	{
		pl.AddTail(ipl[i]);
	}
	return true;
}

int iProtein::GetPL(int num,iPeptide *pl)
{
	if(pl==NULL||num<=0) return 0;
	int tnum=num>count?count:num;
	for(int i=0;i<tnum;i++)
	{
		pl[i]=ipl[i];
	}
	return tnum;
}

bool iProtein::SetAt(int idx,iPeptide &r)
{
	if(idx<0||idx>=count) return false;
	tscore-=ipl[idx].score;
	tscore+=r.score;
	ipl[idx]=r;	
	return true;
}

iPeptide iProtein::GetAt(int idx)
{
   iPeptide itemp;
   if(idx>=0&&idx<count) 
   {
	   itemp=ipl[idx];
   }
   return itemp;
}

bool iProtein::GetAt(int idx,iPeptide &r)
{
	if(idx<0||idx>=count) return false;
	r=ipl[idx];
	return true;
}

iPeptide iProtein::operator [](int idx)
{
   iPeptide itemp;
   if(idx>=0&&idx<count) 
   {
	   itemp=ipl[idx];
   }
   return itemp;
}

iProtein &iProtein::operator =(iProtein &r)
{
	if(count>0&&ipl!=NULL)
	{
		delete []ipl;
		count=0;
	}
	count=r.GetCount();
	if(count>0)
	{
		ipl=new iPeptide[count];
		r.GetPL(count,ipl);
		pAcc=r.pAcc;
		tscore=r.tscore;
		pcount=r.pcount;
		mass=r.mass;
		ratio=r.ratio;
		std=r.std;
	}
	else
	{
		ipl=NULL;
		count=0;
	}
	return *this;
}

int iProtein::GetCount()
{
	return count;
}

void iProtein::UpdateScore()
{
	tscore=0;
	pcount=count;
	for(int i=0;i<count;i++)
	{
		tscore+=ipl[i].score;
	}
}


