// OutRecord.cpp: implementation of the OutRecord class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "OutRecord.h"
#include "afxtempl.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

////////////////////Filter
SFilter::SFilter()
{
	XCorr[0]=1.5f;
	XCorr[1]=2.0f;
	XCorr[2]=2.5f;
	detCn=0.1f;
	RSp=4;
}

SFilter::~SFilter()
{

}

SFilter::SFilter(const SFilter &r)
{
	XCorr[0]=r.XCorr[0];
	XCorr[1]=r.XCorr[1];
	XCorr[2]=r.XCorr[2];
	detCn=r.detCn;
	RSp=r.RSp;
}

SFilter &SFilter::operator=(const SFilter &r)
{
	XCorr[0]=r.XCorr[0];
	XCorr[1]=r.XCorr[1];
	XCorr[2]=r.XCorr[2];
	detCn=r.detCn;
	RSp=r.RSp;
	return *this;
}

bool SFilter::IsPassFilter(SRecord *pep,int ch)
{
	if(pep->detCn<detCn) return false;
	if(ch<=0) return false;
	if(ch>3) ch=3;
	if(pep->XCorr<XCorr[ch-1]) return false;
	if(pep->RankSp>RSp) return false;
	return true;
}

bool SFilter::IsValidate()
{
	if(XCorr[0]<0||XCorr[0]>10) return false;
	if(XCorr[1]<0||XCorr[1]>10) return false;
	if(XCorr[2]<0||XCorr[2]>10) return false;
	if(detCn<0||detCn>0.99999) return false;
	if(RSp<1||RSp>500) return false;
	return true;
}

void SFilter::WriteToFile(CStdioFile &dFile)
{
	CString tmp;
	tmp.Format("XCorr_CH1=%lf\n",XCorr[0]);
	dFile.WriteString(tmp);

	tmp.Format("XCorr_CH2=%lf\n",XCorr[1]);
	dFile.WriteString(tmp);

	tmp.Format("XCorr_CH3=%lf\n",XCorr[2]);
	dFile.WriteString(tmp);

	tmp.Format("detCn=%lf\n",detCn);
	dFile.WriteString(tmp);

	tmp.Format("RSp=%d\n",RSp);
	dFile.WriteString(tmp);
}

bool SFilter::LoadFromFile(CStdioFile &dFile)
{
	CString tmp;
	if(!dFile.ReadString(tmp)) return false;
	int RD=sscanf((LPCTSTR)tmp,"XCorr_CH1=%lf",XCorr);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"XCorr_CH2=%lf",XCorr+1);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"XCorr_CH3=%lf",XCorr+2);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"detCn=%lf",&detCn);
	if(RD!=1) return false;

	if(!dFile.ReadString(tmp)) return false;
	RD=sscanf((LPCTSTR)tmp,"RSp=%d",&RSp);
	if(RD!=1) return false;
	return true;
}
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

OutRecord::OutRecord()
{
	SR=NULL;
	count=0;
	MatchPep=0;
	TotalInt=0;
	EMH=0;
	LowSp=0;
	Charge=0;
}

OutRecord::~OutRecord()
{
	if(count>0&&SR!=NULL)
	{
		delete []SR;
	}
}

OutRecord& OutRecord::operator=(OutRecord &r)
{
	if(count>0&&SR!=NULL)
	{
		delete []SR;
		SR=NULL;
		count=0;
	}
	MatchPep=r.MatchPep;
	TotalInt=r.TotalInt;
	EMH=r.EMH;
	LowSp=r.LowSp;
	Charge=r.Charge;
	count=r.GetCount();
	if(count>0)
	{
		SR=new SRecord[count];
		for(int i=0;i<count;i++)
		{
			r.GetAt(i,&SR[i]);
		}
	}
	return *this;
}

bool OutRecord::RemTailSP(char *str)
{
    int len;
	if(str==NULL) return false;
	len=strlen(str)-1;
	while(len>=0)
	{
		if(str[len]==' ') len--;
		else if(str[len]=='\t') len--;
		else if(str[len]=='\n') len--;
		else if(str[len]=='\r') len--;
		else break;
	}
	str[len+1]='\0';
	return true;
}

bool OutRecord::StrTailFind(char *str,char *substr)
{
	int len;
	char *temp;
	RemTailSP(str);
	len=strlen(str);
	temp=strstr(str,substr);
	if(temp==NULL) return false;
	if(strcmp(temp,substr)==0) return true;
	return false;
}

bool OutRecord::ReadFromFile(FILE *fp,char *outtitle)
{
  bool IsTruboSequest=false;
  char buf[MAX_LINEW];  
  char buf1[MAX_LINEW];
  int n; int i=0;int index,m,k;
  if(fp==NULL) return false; 
  if(feof(fp)) return false;
  sFName=outtitle;
  SRecord rec;
  CList<SRecord,SRecord&> RList;
  while(fgets(buf,MAX_LINEW-1,fp))
  {
	  if(strstr(buf,"(M+H)+ mass =")!=NULL)
	  {
		  n=strlen(buf);
		  index=0;
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  if(index<n)
		  {
             sscanf(buf+index,"%f",&EMH);		
		  }
		  while(index<n&&buf[index]!='+')index++;
		  index++;
		  //Get charge
		  if(index<n)
		  {
			  Charge=buf[index]-'0';
		  }
	  }
	  else if(strstr(buf,"total inten = ")!=NULL)
	  {
		  n=strlen(buf);
		  index=0;
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  i=0;
		  while(index<n&&buf[index]!=',')
		  {
			  buf1[i]=buf[index];
			  i++;index++;
		  }
		  buf1[i]='\0';
		  sscanf(buf1,"%f",&TotalInt);
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  i=0;
		  while(index<n&&buf[index]!=',')
		  {
			  buf1[i]=buf[index];
			  i++;index++;
		  }
		  buf1[i]='\0';
		  sscanf(buf1,"%f",&LowSp);
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  sscanf(buf+index,"%d",&MatchPep);
	  }
	  else if(strlen(buf)<=2) break;
  }
  fgets(buf,MAX_LINEW-1,fp); /*skip the file head*/
  
  if(strstr(buf,"Id")!=NULL) IsTruboSequest=true;
  if(strstr(buf,"deltCn")==NULL) return false;
  if(fgets(buf,MAX_LINEW-1,fp)==NULL) return false;  
  if(strstr(buf,"--------")==NULL) return false;  

  while(!feof(fp))
  { 
	fgets(buf,MAX_LINEW-1,fp);	/*read a line*/
	n=strlen(buf);
	if(n<=4) break;
    index=0;
	while(index<n&&buf[index]==' ')index++;

	if(index<=5)
	{
         for(m=0;m<n; m++)
		 {
            if (buf[m] == '/') buf[m] = ' ';      
		 }	      
		 while(index<n&&buf[index]!=' ') index++;
         while(index<n&&buf[index]==' ') index++;
		 while(index<n&&buf[index]!=' ') index++;
		 while(index<n&&buf[index]==' ') index++;
		 m=0;
		 while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
		 buf1[m]='\0';
         sscanf(buf1,"%d",&rec.RankSp);
         while(index<n&&buf[index]==' ') index++;
	     if(IsTruboSequest)
		 {	  
              while(index<n&&buf[index]!=' ') index++;
		 }		 
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
		 buf1[m]='\0';
	     sscanf(buf1,"%f",&rec.MH);	
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
		 buf1[m]='\0';
         sscanf(buf1,"%f",&rec.detCn);
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';
         sscanf(buf1,"%f",&rec.XCorr);
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';
	     sscanf(buf1,"%f",&rec.Sp);
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
		 while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';	
		 sscanf(buf1,"%d",&rec.mIons); 
		 while(index<n&&buf[index]==' ') index++;
		 m=0;
		 while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';	
	     sscanf(buf1,"%d",&rec.tIons);
 	     while(index<n&&buf[index]==' ') index++;
         m=0;
	     while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';
		 rec.Ref=buf1;
	     while(index<n&&buf[index]==' ') index++;
		 if(buf[index]=='+')                         /*skip +*/
		 { 
		   while(index<n&&buf[index]!=' ') index++;
	       while(index<n&&buf[index]==' ') index++;
		 }
		 m=0;
	     while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}	
         buf1[m]='\0';
		 if(m>4)
		 {
			 if(buf1[1]=='.')
			 {
               for(k=0;k<m-2;k++)buf1[k]=buf1[k+2];
			   while(m>0&&buf1[m]!='.')m--;
			   buf1[m]='\0';			   
			 }
		 }
		 rec.Seq=buf1;
		RList.AddTail(rec);	
	}//if(index<=5)
	else
	{
		if(IsTruboSequest)
		{
			while(index<n&&buf[index]!=' ') index++;//skip the global id
			while(index<n&&buf[index]==' ') index++;
		}
		m=0;
		while(index<n&&buf[index]!=' ') 
		{
			if(buf[index]=='|'||buf[index]=='\n'||buf[index]=='\r')break;
			buf1[m]=buf[index];
			index++;m++;
		}
		buf1[m]='\0';	
		rec.Ref+="+";
		rec.Ref+=buf1;	
	}
  } 
  while(!feof(fp))//seek to next scan
  {
	  fgets(buf,MAX_LINEW-1,fp);	/*read a line*/
	  if(StrTailFind(buf,".out")) 
	  {
          strcpy(outtitle,buf);
		  break;
	  }
  }
  if(count>0&&SR!=NULL)
  {
	  delete []SR;
	  SR=NULL;
  }
  count=RList.GetCount();
  if(count>0)
  {
	  SR=new SRecord[count];
	  POSITION pos=RList.GetHeadPosition();
	  for(i=0;i<count;i++)
	  {
		  SR[i]=RList.GetNext(pos);
	  }
  }
  return true;
}

bool OutRecord::ReadFromFile(CString fname)
{
  FILE *fp;
  fp=fopen(fname,"r");
  if(fp==NULL) return false;
  bool IsTruboSequest=false;
  char buf[MAX_LINEW];  
  char buf1[MAX_LINEW];
  int n; int i=0;int index,m,k;
  sFName=fname;
  SRecord rec;
  CList<SRecord,SRecord&> RList;
  fgets(buf,MAX_LINEW-1,fp);
  while(fgets(buf,MAX_LINEW-1,fp))
  {
	  if(strstr(buf,"(M+H)+ mass =")!=NULL)
	  {
		  n=strlen(buf);
		  index=0;
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  if(index<n)
		  {
             sscanf(buf+index,"%f",&EMH);		
		  }
		  while(index<n&&buf[index]!='+')index++;
		  index++;
		  if(index<n)
		  {
			  Charge=buf[index]-'0';
		  }
	  }
	  else if(strstr(buf,"total inten = ")!=NULL)
	  {
		  n=strlen(buf);
		  index=0;
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  i=0;
		  while(index<n&&buf[index]!=',')
		  {
			  buf1[i]=buf[index];
			  i++;index++;
		  }
		  buf1[i]='\0';
		  sscanf(buf1,"%f",&TotalInt);
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  i=0;
		  while(index<n&&buf[index]!=',')
		  {
			  buf1[i]=buf[index];
			  i++;index++;
		  }
		  buf1[i]='\0';
		  sscanf(buf1,"%f",&LowSp);
		  while(index<n&&buf[index]!='=')index++;
		  index++;
		  sscanf(buf+index,"%d",&MatchPep);
	  }
	  else if(strlen(buf)<=2) break;
  }
  fgets(buf,MAX_LINEW-1,fp); /*skip the file head*/
  
  if(strstr(buf,"Id")!=NULL) IsTruboSequest=true;
  if(strstr(buf,"deltCn")==NULL) {fclose(fp); return false;}
  if(fgets(buf,MAX_LINEW-1,fp)==NULL) {fclose(fp);return false;}
  if(strstr(buf,"--------")==NULL) {fclose(fp);return false;}

  while(!feof(fp))
  { 
	fgets(buf,MAX_LINEW-1,fp);	/*read a line*/
	n=strlen(buf);
	if(n<=4) break;
    index=0;
	while(index<n&&buf[index]==' ')index++;

	if(index<=5)
	{
         for(m=0;m<n; m++)
		 {
            if (buf[m] == '/') buf[m] = ' ';      
		 }	      
		 while(index<n&&buf[index]!=' ') index++;
         while(index<n&&buf[index]==' ') index++;
		 while(index<n&&buf[index]!=' ') index++;
		 while(index<n&&buf[index]==' ') index++;
		 m=0;
		 while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
		 buf1[m]='\0';
         sscanf(buf1,"%d",&rec.RankSp);
         while(index<n&&buf[index]==' ') index++;
	     if(IsTruboSequest)
		 {	  
              while(index<n&&buf[index]!=' ') index++;
		 }		 
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
		 buf1[m]='\0';
	     sscanf(buf1,"%f",&rec.MH);	
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
		 buf1[m]='\0';
         sscanf(buf1,"%f",&rec.detCn);
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';
         sscanf(buf1,"%f",&rec.XCorr);
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
         while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';
	     sscanf(buf1,"%f",&rec.Sp);
	     while(index<n&&buf[index]==' ') index++;
		 m=0;
		 while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';	
		 sscanf(buf1,"%d",&rec.mIons); 
		 while(index<n&&buf[index]==' ') index++;
		 m=0;
		 while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';	
	     sscanf(buf1,"%d",&rec.tIons);
 	     while(index<n&&buf[index]==' ') index++;
         m=0;
	     while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}
         buf1[m]='\0';
		 rec.Ref=buf1;
	     while(index<n&&buf[index]==' ') index++;
		 if(buf[index]=='+')                         /*skip +*/
		 { 
		   while(index<n&&buf[index]!=' ') index++;
	       while(index<n&&buf[index]==' ') index++;
		 }
		 m=0;
	     while(index<n&&buf[index]!=' ') {buf1[m]=buf[index];index++;m++;}	
         buf1[m]='\0';
		 if(m>4)
		 {
			 if(buf1[1]=='.')
			 {
               for(k=0;k<m-2;k++)buf1[k]=buf1[k+2];
			   while(m>0&&buf1[m]!='.')m--;
			   buf1[m]='\0';			   
			 }
		 }
		 rec.Seq=buf1;
		RList.AddTail(rec);	
	}//if(index<=5)
	else
	{
		if(IsTruboSequest)
		{
			while(index<n&&buf[index]!=' ') index++;//skip the global id
			while(index<n&&buf[index]==' ') index++;
		}
		m=0;
		while(index<n&&buf[index]!=' ') 
		{
			if(buf[index]=='|'||buf[index]=='\n'||buf[index]=='\r')break;
			buf1[m]=buf[index];
			index++;m++;
		}
		buf1[m]='\0';	
		rec.Ref+="+";
		rec.Ref+=buf1;	
	}
  }  
  fclose(fp);
  if(count>0&&SR!=NULL)
  {
	  delete []SR;
	  SR=NULL;
  }
  count=RList.GetCount();
  if(count>0)
  {
	  SR=new SRecord[count];
	  POSITION pos=RList.GetHeadPosition();
	  for(i=0;i<count;i++)
	  {
		  SR[i]=RList.GetNext(pos);
	  }
  } 
  return true;
}

int OutRecord::GetCount()
{
	return count;
}

bool OutRecord::GetAt(int idx,SRecord *rec)
{
	if(idx<0||idx>=count) return false;
	*rec=SR[idx];
	return true;
}

bool OutRecord::GetAt(int idx,SRecord &rec)
{
	if(idx<0||idx>=count) return false;
	rec=SR[idx];
	return true;
}

bool OutRecord::SetAt(int idx,SRecord *rec)
{
	if(idx<0||idx>=count) return false;
	SR[idx]=*rec;
	return true;
}

float OutRecord::GetdetCn()
{
	if(count>=2) return SR[1].detCn;
	else return 0;
}
