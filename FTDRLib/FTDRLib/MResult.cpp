// MResult.cpp: implementation of the MResult class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MResult.h"
#include "math.h"
#define _MATRIX_USE_STATIC_LIB 
#include "MSparser\msparser.hpp"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace matrix_science;
////////////////////////Filter
Filter::Filter()
{
	MinPepScore=0.0f;
	DefaultMin=0.0f;
	MaxRank=10;
	MaxMiss=2;
	MET=1.5e-5f;
	IsUseDefault=TRUE;
}

Filter::~Filter()
{
}

Filter::Filter(const Filter &r)
{
	MinPepScore=r.MinPepScore;
	DefaultMin=r.DefaultMin;
	MaxRank=r.MaxRank;
	MaxMiss=r.MaxMiss;
	MET=r.MET;
	IsUseDefault=r.IsUseDefault;	
}

Filter &Filter::operator=(const Filter &r)
{
	MinPepScore=r.MinPepScore;
	DefaultMin=r.DefaultMin;
	MaxRank=r.MaxRank;
	MaxMiss=r.MaxMiss;
	MET=r.MET;
	IsUseDefault=r.IsUseDefault;
	return *this;
}

bool Filter::IsPassFilter(iPeptide *ptemp)
{
	if(ptemp->Sequence.IsEmpty())return false;
	if(IsUseDefault)
	{
		if(ptemp->score<DefaultMin)return false;
	}
	else if(ptemp->score<MinPepScore) return false;
	if(fabs(ptemp->dm)>=ptemp->CalcMass*MET) return false;
	if(ptemp->miss>MaxMiss)  return false;
	if(ptemp->rank>MaxRank)  return false;
	return true;
}

bool Filter::IsValidate()
{
	if(MinPepScore<0) return false;
	if(DefaultMin<0) return false;
	if(MaxRank<1||MaxRank>10) return false;
	if(MaxMiss<0||MaxMiss>10) return false;
	if(MET<0||MET>100) return false;
	return true;
}

void Filter::WriteToFile(CStdioFile &dFile)
{
	CString tmp;
	tmp.Format("MinPepScore=%lf\n",MinPepScore);
	dFile.WriteString(tmp);

	tmp.Format("DefaultMin=%lf\n",DefaultMin);
	dFile.WriteString(tmp);

	tmp.Format("MaxRank=%d\n",MaxRank);
	dFile.WriteString(tmp);

	tmp.Format("MaxMiss=%d\n",MaxMiss);
	dFile.WriteString(tmp);

	tmp.Format("MET=%lf\n",MET);
	dFile.WriteString(tmp);

	tmp.Format("IsUseDefault=%d\n",IsUseDefault);
	dFile.WriteString(tmp);
}

bool Filter::LoadFromFile(CStdioFile &dFile)
{
	CString tmp;
	if(!dFile.ReadString(tmp)) return false;
	int RD=sscanf((LPCTSTR)tmp,"MinPepScore=%lf",&MinPepScore);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"DefaultMin=%lf",&DefaultMin);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"MaxRank=%d",&MaxRank);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"MaxMiss=%d",&MaxMiss);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"MET=%lf",&MET);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"IsUseDefault=%d",&IsUseDefault);
	if(RD!=1) return false;
	return true;
}
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MResult::MResult()
{
	count=0;
	PEPLIST=NULL;	
}

MResult::~MResult()
{
	if(PEPLIST!=NULL&&count>0) delete []PEPLIST;

}

int MResult::GetCount()
{
	return count;
}

MResult::MResult(int num)
{
	if(num>0)
	{
		PEPLIST=new iPeptide[num];
		count=num;
	}
	else
	{
		PEPLIST=NULL;
		count=0;
	}
}

bool MResult::SetPL(int num,iPeptide *pl)
{
	if(PEPLIST!=NULL&&count>0) delete []PEPLIST;
	if(num<=0||pl==NULL)
	{
		PEPLIST=NULL;
		count=0;
		return false;
	}
	PEPLIST=new iPeptide[num];
	count=num;
	for(int i=0;i<num;i++)PEPLIST[i]=pl[i];
	return true;
}

bool MResult::SetPL(CList<iPeptide,iPeptide&> *pl)
{
	if(PEPLIST!=NULL&&count>0) 
	{
		delete []PEPLIST;
		count=0;	
		PEPLIST=NULL;
	}
	if(pl==NULL)return false;
	int num=pl->GetCount();
	if(num<=0)	return false;
	POSITION pos=pl->GetHeadPosition();
	PEPLIST=new iPeptide[num];			
	count=num;
	for(int i=0;i<num;i++)PEPLIST[i]=pl->GetNext(pos);
	return true;
}

bool MResult::SetPL(CList<iPeptide,iPeptide&> &pl)
{
	if(PEPLIST!=NULL&&count>0) delete []PEPLIST;
	int num=pl.GetCount();
	if(num<=0)
	{
		PEPLIST=NULL;
		count=0;
		return false;
	}
	POSITION pos=pl.GetHeadPosition();
	PEPLIST=new iPeptide[num];			
	count=num;
	for(int i=0;i<num;i++)
		PEPLIST[i]=pl.GetNext(pos);
	return true;
}

bool MResult::GetPL(CList<iPeptide,iPeptide&> *pl)
{
	if(pl==NULL) return false;
	pl->RemoveAll();
	for(int i=0;i<count;i++)pl->AddTail(PEPLIST[i]);
	return true;
}

bool MResult::GetPL(CList<iPeptide,iPeptide&> &pl)
{
	pl.RemoveAll();
	for(int i=0;i<count;i++)pl.AddTail(PEPLIST[i]);
	return true;
}

int MResult::GetPL(int num,iPeptide *pl)
{
	if(pl==NULL||num<=0) return 0;
	int tnum=num>count?count:num;
	for(int i=0;i<tnum;i++)pl[i]=PEPLIST[i];
	return tnum;
}

bool MResult::SetAt(int idx,iPeptide &r)
{
	if(idx>=count||idx<0) return false;
	PEPLIST[idx]=r;
	return true;
}

bool MResult::SetAt(int idx,iPeptide *r)
{
	if(idx>=count||idx<0) return false;
	PEPLIST[idx]=*r;
	return true;
}

bool MResult::GetAt(int idx,iPeptide &r)
{
	if(idx>=count||idx<0) return false;
	r=PEPLIST[idx];
	return true;
}

bool MResult::GetAt(int idx,iPeptide *r)
{
	if(idx>=count||idx<0) return false;
	*r=PEPLIST[idx];
	return true;
}

MResult &MResult::operator =(MResult &r)
{
	if(PEPLIST!=NULL&&count>0) delete []PEPLIST;
	int num=r.GetCount();
	if(num<=0)
	{
		PEPLIST=NULL;
		count=0;	
	}
	else
	{	
		PEPLIST=new iPeptide[num];			
		count=num;
		r.GetPL(count,PEPLIST);	
	}
	cFilter=r.cFilter;
	return *this;
}

bool MResult::ReadOnePep(CString sTemp,iPeptide &ptemp)
{
	ptemp.Clear();
	CString item;
	int sidx=sTemp.Find("</A>");
	if(sidx!=-1)
	{
	   sidx+=4;
	   if(sidx>=sTemp.GetLength()-1) return false;
	   sTemp=sTemp.Mid(sidx);	  	
	}
	else return false;
	sTemp.Replace("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;","\t");
	sTemp.Replace("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;","\t");
	sTemp.Replace("&nbsp;&nbsp;&nbsp;&nbsp;","\t");
	sTemp.Replace("&nbsp;&nbsp;&nbsp;","\t");
	sTemp.Replace("&nbsp;&nbsp;","\t");
	sTemp.Replace("&nbsp;","\t");
	sTemp.Replace("\t\t;","\t");

	//sidx=sTemp.Find("&nbsp;");
	//while(sidx!=-1)
	//{
	//	sTemp.Delete(sidx,6);
	//	sidx=sTemp.Find("&nbsp;");
	//}
	sidx=sTemp.Find("<I>");
	if(sidx!=-1) sTemp.Delete(sidx,3);
	sidx=sTemp.Find("<B>");
	if(sidx!=-1) 
	{
		ptemp.IsBold=true;
		sTemp.Delete(sidx,3);
	}
	sidx=sTemp.Find("<FONT       color=#ff0000>");
	if(sidx==-1)
	{
		sidx=sTemp.Find("<FONT       COLOR=#ff0000>");
	}
	if(sidx!=-1)
	{
		sidx+=26;
		ptemp.color=COLOR_RED;
		if(sidx>=sTemp.GetLength()-1) return false;
		sTemp=sTemp.Mid(sidx);
	}
	else 
	{
		sidx=sTemp.Find('\t');
		if(sidx>=sTemp.GetLength()-1) return false;
		sTemp=sTemp.Mid(sidx+1);//skip the index
		ptemp.color=COLOR_BLACK;
   }
	//cout<<"New line:"<<(LPCTSTR)sTemp<<endl;
	sidx=sTemp.Find('\t');
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);
	sscanf(item,"%f",&ptemp.Observed);
	//ptemp.Observed=atof(item);
	if(sidx>=sTemp.GetLength()-1) return false;
	sTemp=sTemp.Mid(sidx+1);

	sidx=sTemp.Find('\t');
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);
	sscanf(item,"%f",&ptemp.CalcMass);
	//ptemp.CalcMass=atof(item);
	if(sidx>=sTemp.GetLength()-1) return false;
	sTemp=sTemp.Mid(sidx+1);

	sidx=sTemp.Find('\t');
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);
	sscanf(item,"%f",&ptemp.ObsrvMass);
	//ptemp.CalcMass=atof(item);
	if(sidx>=sTemp.GetLength()-1) return false;
	sTemp=sTemp.Mid(sidx+1);

	sidx=sTemp.Find('\t');//dm
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);
	sscanf(item,"%f",&ptemp.dm);
	//ptemp.da=atof(item);
	if(sidx>=sTemp.GetLength()-1) return false;
	sTemp=sTemp.Mid(sidx+1);

	sTemp.Remove('(');
	sTemp.Remove(')');
	sidx=sTemp.Find('\t');//miss
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);
	sscanf(item,"%d",&ptemp.miss);	
	//ptemp.miss=atoi(item);
	if(sidx>=sTemp.GetLength()-1) return false;
	sTemp=sTemp.Mid(sidx+1);

	sidx=sTemp.Find('\t');//score
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);
	sscanf(item,"%d",&ptemp.score);
	//ptemp.score=atoi(item);
	if(sidx>=sTemp.GetLength()-1) return false;
	sTemp=sTemp.Mid(sidx+1);

	sidx=sTemp.Find('\t');//rank
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);
	sscanf(item,"%d",&ptemp.rank);
	if(sidx>=sTemp.GetLength()-1) return false;
	sTemp=sTemp.Mid(sidx+1);

	sidx=sTemp.Find("</");
	if(sidx==-1) return false;
	item=sTemp.Mid(0,sidx);	
	sidx=item.Find('+');
	if(sidx==-1) ptemp.Sequence=item;
	else
	{
		ptemp.Sequence=item.Mid(0,sidx);
		ptemp.Modif=item.Mid(sidx+1);
	}
	return true;
}

int MResult::ReadFromFile(LPCTSTR fname)
{
	CStdioFile dFile;
	if(!dFile.Open(fname,CFile::modeRead|CFile::typeText))
	{
			//AfxMessageBox("Can not open this file!");
			return 0;
	}
	IsCarb=false;
	CString sTemp;
	bool IsLoad=false;
	bool IsReadPep=false;
	CList<iPeptide,iPeptide&> pepl;	
	iPeptide ptemp;
	CString item;
	if(cFilter.IsUseDefault)
	{
		while(dFile.ReadString(sTemp))
		{
			if(sTemp.Find("<H3>Probability Based Mowse Score</H3>")!=-1)
			{
				item+=sTemp;
				break;
			}

			//if(sTemp.Find("-10*Log(P)")!=-1) break;
		}
		while(item.Find("indicate identity or extensive homology")==-1)
		{	
			if(!dFile.ReadString(sTemp)) 
			{
				dFile.Close();
				return 0;
			}
			item+=sTemp;
		}

		int sidx=item.Find("&gt; ");
		if(sidx!=-1)
		{
			sidx+=5;
			item=item.Mid(sidx);
			sidx=item.Find(" ");
			if(sidx>0)	item=item.Mid(0,sidx);
			int tempi;
			sscanf(item,"%d",&tempi);
			cFilter.DefaultMin=(float)tempi;
		}
	}
	while(dFile.ReadString(sTemp))
	{
		CString str=sTemp;
		str.MakeLower();	
		if(sTemp.Find("<TD><INPUT type=checkbox CHECKED")!=-1)
		{
			ptemp.IsChecked=true;
		}
		else if(sTemp.Find("value=\"")!=-1)
		{
			int tidx=sTemp.Find("title(");
			if(tidx!=-1)
			{
				tidx+=6;
				if(tidx<sTemp.GetLength())
				{
					CString item=sTemp.Mid(tidx);
					tidx=item.Find(")");
					if(tidx!=-1)
					{
						item=item.Mid(0,tidx);
						item.Replace("%2e",".");
						ptemp.ScanNumber=GetScanNumber(item);
						ptemp.fName=item;
					}
				}
			}
		}
		else if(str.Find("target=_blank>")!=-1)
		{
			CString stemp1=sTemp;
			while(sTemp.Find("</TR>")==-1)
			{
				if(!dFile.ReadString(sTemp)) break;
				if(sTemp.Find("<TR>")!=-1) break;
				stemp1+=sTemp;
			}			
			if(ReadOnePep(stemp1,ptemp))
			{
				if(cFilter.IsPassFilter(&ptemp))			
				{						
					pepl.AddTail(ptemp);
					count++;
				}			
			}
			ptemp.IsChecked=false;
			ptemp.ScanNumber=0;
		}
		else if(sTemp.Find("Unassigned queries:")!=-1) break;
	}
	while(dFile.ReadString(sTemp))
	{
		if(sTemp.Find("<H3>Search Parameters</H3>")!=-1)
		{
			while(dFile.ReadString(sTemp))
			{
				if(sTemp.Find("<B>Fixed modifications ")!=-1)
				{
					if(sTemp.Find("Carbamidomethyl (C)")!=-1)
						IsCarb=true;				
					break;
				}
			}
		}
	}	
	dFile.Close();
	RemoveRon(&pepl);
	SetPL(&pepl);
	return count;
}

bool MResult::RemoveRon(CList<iPeptide,iPeptide&> *ipl)
{
	if(ipl==NULL) return false;
	int count=ipl->GetCount();
	if(count<=0) return false;
	iPeptide *pltemp;
	iPeptide itemp;
	pltemp=new iPeptide[count];
	POSITION pos=ipl->GetHeadPosition();
	pltemp[0]=ipl->GetNext(pos);
	int i,k,num=1;	
	for(i=1;i<count;i++)
	{
		itemp=ipl->GetNext(pos);
	    int ch2=(int)((itemp.ObsrvMass+1.007825)/itemp.Observed+0.3);
		for(k=0;k<num;k++)
		{
			if(pltemp[k].fName==itemp.fName) break;
		}
		if(k==num)
		{
			pltemp[num]=itemp;
			num++;
		}
	}
	ipl->RemoveAll();
	for(i=0;i<num;i++)ipl->AddTail(pltemp[i]);
	delete []pltemp;
	return true;
}

void MResult::Clear()
{
	if(PEPLIST!=NULL&&count>0) delete []PEPLIST;
	count=0;
	PEPLIST=NULL;
}

int MResult::GetScanNumber(CString item)
{
	int idx=item.Find("FinneganScanNumber");
	if(idx!=-1)
	{
		idx+=18;
		item=item.Mid(idx);
		int tmpi;
		sscanf((LPCTSTR)item,"%d",&tmpi);
		return tmpi;
	}
	int len=item.GetLength();
	if(len<3) return 0;
	CString tstr;
	tstr=item.Mid(len-3);
	if(tstr=="out"||tstr=="dta") 
	{
		len-=4;
		item=item.Mid(0,len);
	}
	item=item.Mid(0,len-2);//remove charge	
	int i=len-3;
	while(i>=0&&item[i]!='.') i--;
	if(i<0) return 0;
	i++;
	item=item.Mid(i,len-2-i);
	return atoi(item);
}


//	3232	702.4007	1402.7869	1402.6953	0.0916	1	3	1.9e+002	3	R.IFAENNTARDPR.L
bool MResult::ReadOnePep2(CString sTemp,iPeptide &ptemp)
{
	ptemp.Clear();
	CString str=sTemp;
	str.MakeLower();
	int idx=str.Find("onmouseout=");
	if(idx==-1) return false;
	if(str.Find("color=#ff0000",idx)!=-1) ptemp.color=COLOR_RED;
	if(sTemp.Find("<B>",idx)!=-1) ptemp.IsBold=true;
	RemoveMarker(sTemp);
	sTemp.Remove(' ');
	sTemp.Remove('\t');
	ReplaceTokens(sTemp);
	sTemp.Remove('(');
	sTemp.Remove(')');
	int QId;
	char seqbuf[1024];
	int RDNum=sscanf(sTemp,"\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%e\t%d\t%s",\
		&QId,&ptemp.Observed,&ptemp.ObsrvMass,&ptemp.CalcMass,&ptemp.dm,&ptemp.miss,\
		&ptemp.score,&ptemp.Expect,&ptemp.rank,seqbuf);	
	if(Isppm) ptemp.dm=ptemp.dm*ptemp.CalcMass/1e6;
	if(RDNum!=10) return false;
	ptemp.Sequence=seqbuf;
	return true;
}

void MResult::ReplaceTokens(CString &sTemp)
{
	int slen=sTemp.GetLength();
	int len=sTemp.GetLength();
	CString str;
	bool IsMarker=false;
	bool IsAdd=false;
	for(int i=0;i<len;i++)
	{
		if(sTemp[i]=='&')
		{
			IsMarker=true;
			if(!IsAdd)
			{
				str+="\t";
				IsAdd=true;
			}
		}
		else if(sTemp[i]==';')IsMarker=false;
		else if(!IsMarker)
		{
			str+=sTemp[i];	
			if(IsAdd) IsAdd=false;
		}
	}
	sTemp=str;
}

void MResult::RemoveSPs(CString &sTemp)
{
	sTemp.Replace("%2e",".");
	int len=sTemp.GetLength();
	CString str;
	bool IsMarker=false;
	int i=0;
	while(i<len)
	{
		if(sTemp[i]=='%')
		{
			if(!IsMarker)sTemp+=' ';			
			IsMarker=true;
			i+=3;//skip the three
			if(i>=len) break;
		}
		else
		{
			str+=sTemp[i];
			i++;
			IsMarker=false;
		}
	}
	sTemp=str;
}

int MResult::SubStrCount(CString &str,char *substr)
{
	char *pStr;
	pStr=str.GetBuffer();
	int slen=str.GetLength();
	int count=0;
    pStr=strstr(pStr,substr);
	while(pStr!=NULL)
	{
		count++;
		pStr=strstr(pStr,substr);
	}
	str.ReleaseBuffer();
	return count;
}

void MResult::MakeLower(CString &str)
{
	int slen=str.GetLength();
	CString Stemp;
	for(int i=0;i<slen;i++)
	{
		char c=str.GetAt(i);
		if(c>='A'&&c<='Z') c=c-'A'+'a';
		Stemp+=c;
	}
	str=Stemp;
}

bool MResult::ReadAtable(CStdioFile &dFile,CString &str)
{
	str.Empty();
	CString tmpStr;	
	bool RDS=false;

	while(dFile.ReadString(tmpStr))
	{
		if(tmpStr.Find("<TR>")!=-1) 
		{
			RDS=true;
			break;
		}
	}
	if(!RDS) return false;
	str=tmpStr;
	int idx=str.Find("</TR>");
	if(idx!=-1)
	{
		str=str.Mid(0,idx+4);
		return true;
	}
	while(dFile.ReadString(tmpStr))
	{
		idx=tmpStr.Find("</TR>");
		if(idx!=-1) 
		{
			str+=tmpStr.Mid(0,idx+4);
			return true;			
		}
		str+=tmpStr;
	}
	return false;	
}

int MResult::ReadFromFile2(LPCTSTR fname)
{
	CStdioFile dFile;
	IsCarb=false;
	if(!dFile.Open(fname,CFile::modeRead|CFile::typeText))
	{
			//AfxMessageBox("Can not open this file!");
			return 0;
	}
	Isppm=false;
	CString sTemp,tmpStr;
	int total_num=0;
	while(ReadAtable(dFile,sTemp))
	{
		sTemp.MakeLower();
		if(sTemp.Find("<a name=hit")!=-1||sTemp.Find("<a name=\"hit")!=-1)
		{
			total_num = total_num+1;
			if(total_num>=2) break;
		}
		if(sTemp.Find("observed")!=-1&&sTemp.Find("ppm")!=-1) 
		{
			Isppm=true;
			break;
		}	
	}

	dFile.SeekToBegin();

	CList<iPeptide,iPeptide&> pepl;	 
	iPeptide ptemp;
	if(cFilter.IsUseDefault)
	{
		while(dFile.ReadString(sTemp))
		{
			if(sTemp.Find("-10*Log(P)")!=-1) break;
		}
		tmpStr=sTemp;
		while(sTemp.Find("<P>")==-1)
		{
			dFile.ReadString(sTemp);
			tmpStr+=sTemp;
		}	
		int sidx=tmpStr.Find("&gt; ");
		if(sidx!=-1)
		{
			sidx+=5;
			tmpStr=tmpStr.Mid(sidx);
			sidx=tmpStr.Find(" ");
			if(sidx>0)	tmpStr=tmpStr.Mid(0,sidx);
			int tempi;
			sscanf(tmpStr,"%d",&tempi);
			cFilter.DefaultMin=(float)tempi;
		}
	}

	CString str;
	while(ReadAtable(dFile,str))
	{				;
		str.MakeLower();
		if(str.Find("<a name=hit")!=-1||str.Find("<a name=\"hit")!=-1) break; //seek to the  protein record line		
	}	

	while(ReadAtable(dFile,sTemp))
	{	
		str=sTemp;
		str.MakeLower();				
		if((str.Find("type=checkbox")!=-1||str.Find("type=\"checkbox")!=-1)&&(str.Find("checked")!=-1))			
		{	
			ptemp.IsChecked=true;
			int tidx=str.Find(" title(");			
			if(tidx!=-1)
			{
				tidx+=5;
				if(tidx<str.GetLength())
				{
					CString item=str.Mid(tidx);
					tidx=item.Find(")");
					if(tidx!=-1)
					{
						item=item.Mid(0,tidx);
					    item.Replace("%2e",".");						
						RemoveSPs(item);
						ptemp.ScanNumber=GetScanNumber(item);
						ptemp.fName=item;
					}
				}
			}
			if(ReadOnePep2(sTemp,ptemp))
			{
				if(cFilter.IsPassFilter(&ptemp))
				{
					pepl.AddTail(ptemp);
					count++;
					//printf("%d\n",count);
				}								
			}
			ptemp.IsChecked=false;
			ptemp.ScanNumber=0;
			ptemp.fName.Empty();
		}
		else if(sTemp.Find("Peptide matches not assigned to protein hits:")!=-1) break;
		//str.ReleaseBuffer();
	}
	while(dFile.ReadString(sTemp))
	{
		if(sTemp.Find("<H3>Search Parameters</H3>")!=-1)
		{
			while(dFile.ReadString(sTemp))
			{
				if(sTemp.Find("<B>Fixed modifications ")!=-1)
				{
					if(sTemp.Find("Carbamidomethyl (C)")!=-1)
						IsCarb=true;				
					break;
				}
			}
		}
	}	
	dFile.Close();	
	RemoveRon(&pepl);
	SetPL(&pepl);
	return count;
}


int MResult::ReadFromFileNew(LPCTSTR fname)
{
	CStdioFile dFile;
	IsCarb=false;
	if(!dFile.Open(fname,CFile::modeRead|CFile::typeText))
	{
			//AfxMessageBox("Can not open this file!");
			return 0;
	}
	CString sTemp,item;
	bool IsLoad=false;
	bool IsReadPep=false;
	CList<iPeptide,iPeptide&> pepl;
	iPeptide ptemp;
	if(cFilter.IsUseDefault)
	{
		while(dFile.ReadString(sTemp))
		{
			if(sTemp.Find("<H3>Probability Based Mowse Score</H3>")!=-1)
			{
				item+=sTemp;
				break;
			}

			//if(sTemp.Find("-10*Log(P)")!=-1) break;
		}
		while(item.Find("indicate identity or extensive homology")==-1)
		{	
			if(!dFile.ReadString(sTemp)) 
			{
				dFile.Close();
				return 0;
			}
			item+=sTemp;
		}

		int sidx=item.Find("&gt; ");
		if(sidx!=-1)
		{
			sidx+=5;
			item=item.Mid(sidx);
			sidx=item.Find(" ");
			if(sidx>0)	item=item.Mid(0,sidx);
			int tempi;
			sscanf(item,"%d",&tempi);
			cFilter.DefaultMin=(float)tempi;
		}
	}	
	while(dFile.ReadString(sTemp))
	{
		CString str=sTemp;
		str.MakeLower();
		if(sTemp.Find("<INPUT CHECKED name=QUE type=checkbox")!=-1)
		{
			ptemp.IsChecked=true;
		}	
		else if(sTemp.Find("value=\"")!=-1)
		{
			int tidx=sTemp.Find("title(");
			if(tidx!=-1)
			{
				tidx+=6;
				if(tidx<sTemp.GetLength())
				{
					CString item=sTemp.Mid(tidx);
					tidx=item.Find(")");
					if(tidx!=-1)
					{
						item=item.Mid(0,tidx);
						item.Replace("%2e",".");
						ptemp.ScanNumber=GetScanNumber(item);
						ptemp.fName=item;
					}
				}
			}
		}
		else if(sTemp.Find("peptide_view.pl?file=")!=-1)
		{
			while(dFile.ReadString(sTemp))
			{
				if(sTemp.Find("</TT></TD>")!=-1) break;
			}		
			CString stemp1;
			ptemp.Sequence.Empty();
			if(dFile.ReadString(sTemp))
			{	
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}				
				if(stemp1.Find("color=#ff0000")!=-1) ptemp.color=COLOR_RED;
				else ptemp.color=COLOR_BLACK;
				if(stemp1.Find("<B>")!=-1) ptemp.IsBold=true;
				else ptemp.IsBold=false;
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.Observed=(float)atof(stemp1);			
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}	
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.ObsrvMass=(float)atof(stemp1);				
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}	
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.CalcMass=(float)atof(stemp1);			
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}	
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.dm=(float)atof(stemp1);
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}	
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.miss=atoi(stemp1);
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}	
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				stemp1.Remove('(');
				stemp1.Remove(')');
				ptemp.score=atoi(stemp1);
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}	
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.Expect=(float)atof(stemp1);
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}	
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.rank=atoi(stemp1);
			}
			if(dFile.ReadString(sTemp))
			{
				stemp1=sTemp;
				while(sTemp.Find("</TT></TD>")==-1)
				{	
					if(!dFile.ReadString(sTemp)) break;
					stemp1+=sTemp;
				}
				RemoveMarker(stemp1);
				stemp1.Replace("&nbsp;"," ");
				stemp1.Remove(' ');
				ptemp.Sequence=stemp1;
			}

			if(cFilter.IsPassFilter(&ptemp))
			{
				pepl.AddTail(ptemp);
				count++;
			}		
			ptemp.IsChecked=false;
			ptemp.ScanNumber=0;				
		}	
		else if(sTemp.Find("Peptide matches not assigned to protein hits:")!=-1) break;
	}
	while(dFile.ReadString(sTemp))
	{
		if(sTemp.Find("<H3>Search Parameters</H3>")!=-1)
		{
			while(dFile.ReadString(sTemp))
			{
				if(sTemp.Find("<B>Fixed modifications ")!=-1)
				{
					if(sTemp.Find("Carbamidomethyl (C)")!=-1)
						IsCarb=true;				
					break;
				}
			}
		}
	}	
	dFile.Close();	
	RemoveRon(&pepl);
	SetPL(&pepl);
	return count;
}

void MResult::RemoveMarker(CString &sTemp)
{
	int len=sTemp.GetLength();
	CString str;
	bool IsMarker=false;
	for(int i=0;i<len;i++)
	{
		if(sTemp[i]=='<')IsMarker=true;
		else if(sTemp[i]=='>')IsMarker=false;
		else if(!IsMarker)str+=sTemp[i];		
	}
	sTemp=str;
}

int MResult::GetVersion(CString fname)
{
	if(ValidateEXT(fname,"dat")) return 3;
	CStdioFile dFile;
	if(!dFile.Open(fname,CFile::modeRead|CFile::typeText)) return -1;
	CString sTemp;
	int type=0;
	while(dFile.ReadString(sTemp))
	{
		
		if(sTemp.Find("bgColor=#eeeeff noWrap>Significance threshold p&lt;")!=-1) 
			type=1;		
		else if(sTemp.Find("noWrap bgColor=#eeeeff>Significance threshold p&lt;")!=-1) 
			type=1;
		else if(sTemp.Find("<B>emPAI:</B>")!=-1) 
		{
			type=2;
			break;
		}		
	}
	dFile.Close();
	return type;
}

bool MResult::LoadData(CString fname)
{
	int type=GetVersion(fname);
	if(type==1) 
	{
		if(ReadFromFileNew(fname)>0)	return true;
	}
	else if(type==0) 
	{
		if(ReadFromFile(fname)>0)	return true;
	}
	else if(type==2)
	{
		if(ReadFromFile2(fname)>0)	return true;
	}
	else if(type==3) 
	{
		if(LoadFromDat(fname)) return true;
	}
	return false;
}

int MResult::LoadFromDat(LPCTSTR fname)
{
	ms_mascotresfile file(fname);
	if(!file.isValid()) return false;
	if(!file.isMSMS()) return false;
	iPeptide iPep;	
	CList<iPeptide,iPeptide&> pepl;	
	ms_mascotresults * results;
	const int flag= ms_mascotresults::MSRES_GROUP_PROTEINS|ms_mascotresults::MSRES_SHOW_SUBSETS;
	results = new ms_peptidesummary(file,flag,0,20,0,0);	
	if(file.getLastError())
	{
		delete results;
		return false;
	}
	if(cFilter.IsUseDefault) cFilter.DefaultMin=results->getAvePeptideIdentityThreshold(20);
	int total_num=results->getNumberOfHits();
	if(total_num<=0)
	{
		delete results;
		return 0;
	}
	int count = 0;
	//char tmpbuf[100];
	//sprintf(tmpbuf,"Open file successful,parser the %d proteins:",total_num);
	for (int query=1; query <= file.getNumQueries(); query++)
	{
		double threshold = results->getPeptideIdentityThreshold(query, 20);
		ms_peptide * pep;
		if(results->getPeptide(query,1,pep)&&(pep->getIonsScore()>=threshold))//ingrone the unsiginificent matches
		{
			iPep.color=COLOR_RED;
			ms_inputquery *quy=new ms_inputquery(file,query);
			CString title=quy->getStringTitle(true).c_str();			
			delete quy;
			iPep.ScanNumber= GetScanNumber(title);
			iPep.fName=title;
			iPep.score=pep->getIonsScore();		
			iPep.CalcMass=pep->getMrCalc();
			iPep.dm=pep->getDelta();
			iPep.Observed=pep->getObserved();
			iPep.ObsrvMass=pep->getMrExperimental();
			iPep.rank=pep->getRank();
			iPep.miss=pep->getMissedCleavages();

			iPep.Expect=results->getPeptideExpectationValue(iPep.score,query);
			iPep.IsChecked=true;
			iPep.IsBold=true;
			iPep.Sequence=pep->getPeptideStr().c_str();
			//iPep.Modif=pep->getVarModsStr().c_str();
			iPep.Modif=results->getReadableVarMods(query,iPep.rank).c_str();
			if(cFilter.IsPassFilter(&iPep))
			{
				pepl.AddTail(iPep);
				count++;
			}
		}
	}	
	delete results;
	RemoveRon(&pepl);
	SetPL(&pepl);
	return count;
}

bool MResult::ValidateEXT(LPCTSTR fname,const char *EXT)
{
	char tmpFName[MAX_PATH];
	strcpy(tmpFName,fname);	
	int slen=strlen(tmpFName)-1;
	while(slen>=0&&tmpFName[slen]!='.')slen--;
	if(slen==-1) return false;
	strlwr(tmpFName+slen);//make lower case
	if(strstr(tmpFName+slen+1,EXT)!=NULL) return true;
	else return false;
}
