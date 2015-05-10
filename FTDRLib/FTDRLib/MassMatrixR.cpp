#include "StdAfx.h"
#include "MassMatrixR.h"

MassMatrixR::MassMatrixR(void)
{
}

MassMatrixR::~MassMatrixR(void)
{
}


bool MassMatrixR::freadline(FILE *fp,char *buf,int maxcount)
{
	if(fgets(buf,maxcount,fp)==NULL) return false;
	int slen=strlen(buf);
	if(slen>=2&&buf[slen-2]=='\r') buf[slen-2]='\0';
	else if(slen>=1&&buf[slen-1]=='\n') buf[slen-1]='\0';
	return true;
}

char *MassMatrixR::seekTable(int n,char *buf)
{
	char *pStr=strstr(buf,",");
	int i=0;
	while(pStr!=NULL)
	{
		i++;
		if(i==n) break;
		pStr++;
		pStr=strstr(pStr,",");
	}
	if(pStr!=NULL) pStr++;
	return pStr;
}

//Index 0
//scan# 1
//charge 2
//score 3
//pp 4
//pp2 5
//pp_tag 6
//RT_conf 7
//RT(obs) 8
//RT(pred) 9
//PA 10
//m/z 11
//MW(obs) 12
//MW 13
//delta 14
//miss 15
//peptide+modif 16
//protein 17
//2867,2464,3,425,58.4,40.2,3.2,N/A,N/A,N/A,N/A,705.63,2114.87,2114.87,0.00,0,
//MVNNGHSFNVEYDDSQDK OxiM(1),"P00921	|CAH2_BOVIN	CARBONIC	ANHYDRASE	II	(EC	4.2.1.1)	(CARBONATE	DEHYDRATASE	II)	(CA-II)	-	Bos	taurus	(Bovine)."
int MassMatrixR::LoadData(char *fname,vector<Pre_calData> &PL)
{
	PL.clear();
	FILE *fp;
	fp=fopen(fname,"r");
	if(fp==NULL) return false;
	char buf[4096];
	Pre_calData pt;

	char *pStr;
	if(!freadline(fp,buf,4095))
	{
		fclose(fp);
		return false;
	}
	int CT=0;
	while(!feof(fp))
	{
		freadline(fp,buf,4095);
		pStr=seekTable(1,buf);
		int RD=sscanf(pStr,"%d",&(pt.scan));
		if(RD!=1) continue;

		pStr=seekTable(2,buf);
		RD=sscanf(pStr,"%d",&(pt.charge));
		if(RD!=1) continue;

		pStr=seekTable(12,buf);
		RD=sscanf(pStr,"%lf",&(pt.mze));
		if(RD!=1) continue;

		pStr=seekTable(13,buf);
		RD=sscanf(pStr,"%lf",&(pt.mzt));
		if(RD!=1) continue;

		pt.mze+=(pt.charge-1)*1.007296;
		pt.mze/=pt.charge;

		pt.mzt+=(pt.charge-1)*1.007296;
		pt.mzt/=pt.charge;
		PL.push_back(pt);
		CT++;
	}
	fclose(fp);
	return CT;
}