#include "StdAfx.h"
#include "PepXML.h"

PepXML::PepXML(void)
{
	conf_level=1;
	//for(int i=0;i<9;i++)
	//{
	//	ErrTable[i]=0;
	//}
	Min_Prob=0.05;
}

PepXML::~PepXML(void)
{

}

void PepXML::RemTail(char *buf)
{
	int len=strlen(buf);
	if(len<2) return;
	if(buf[len-2]=='\r')buf[len-2]='\0';
	else buf[len-1]=0;	
}
/*
<error_point error="0.000" min_prob="0.87" num_corr="79" num_incorr="0"/>
<error_point error="0.010" min_prob="0.87" num_corr="84" num_incorr="1"/>
<error_point error="0.020" min_prob="0.87" num_corr="90" num_incorr="2"/>
<error_point error="0.025" min_prob="0.87" num_corr="94" num_incorr="3"/>
<error_point error="0.030" min_prob="0.87" num_corr="98" num_incorr="3"/>
<error_point error="0.040" min_prob="0.86" num_corr="107" num_incorr="5"/>
<error_point error="0.050" min_prob="0.85" num_corr="118" num_incorr="6"/>
<error_point error="0.075" min_prob="0.81" num_corr="145" num_incorr="12"/>
<error_point error="0.100" min_prob="0.75" num_corr="171" num_incorr="19"/>
*/
int PepXML::LoadData(char *fname,vector<Pre_calData> &PL)
{
	FILE *fp;
	char buf[4096];
	fp=fopen(fname,"r");
	if(fp==NULL) return 0;

	Pre_calData pt;	
	double st;
	int count=0;
	//seek to the error table
	char *pStr;
	while(!feof(fp))
	{
		fgets(buf,4095,fp);
		pStr=strstr(buf,"<error_point error=");
		if(pStr!=NULL) break;
	}
	//int i;
	while(pStr!=NULL)
	{
		pStr+=20;
		double fdr,minp;
		sscanf(pStr,"%lf",&fdr);
		pStr=strstr(buf,"min_prob=");
		pStr+=10;
		sscanf(pStr,"%lf",&minp);
		ErrTable.push_back(minp);
		FDRTable.push_back(fdr);
		fgets(buf,4095,fp);
		pStr=strstr(buf,"<error_point error=");
	}

	GetCutoff();
	
	//Get scan number and others
	while(!feof(fp))
	{
		fgets(buf,4095,fp);
		if(strstr(buf,"<spectrum_query spectrum=")!=NULL) 
		{
			pStr=strstr(buf,"start_scan=");
			pStr+=12;
			sscanf(pStr,"%d",&pt.scan);		
			
			//not so hig accurate
			//pStr=strstr(buf,"precursor_neutral_mass=");
			//pStr+=24;
			//sscanf(pStr,"%lf",&pt.mze);	

			pStr=strstr(buf,"assumed_charge=");
			pStr+=16;
			pt.charge=*pStr-'0';
		}	
		else if(strstr(buf, "<search_hit hit_rank=\"1\"")!=NULL)
		{
			pStr=strstr(buf,"calc_neutral_pep_mass=");
			pStr+=23;
			sscanf(pStr,"%lf",&pt.mzt);
			pStr=strstr(buf,"massdiff=");
			pStr+=10;
			sscanf(pStr,"%lf",&pt.mze);
			pt.mze+=pt.mzt;
		}
		else if(strstr(buf,"<peptideprophet_result probability=")!=NULL)
		{
			pStr=strstr(buf,"probability=");
			pStr+=13;
			sscanf(pStr,"%lf",&st);
			//rt.scores.insert(rt.scores.begin(),st);
		}
		else if(strstr(buf,"</search_score_summary>")!=NULL)
		{
			//pt.RD.push_back(rt);
			if(pt.charge>0)
			{
				pt.mze=(pt.mze+1.00729*pt.charge)/pt.charge;
				pt.mzt=(pt.mzt+1.00729*pt.charge)/pt.charge;
				if(st>Min_Prob)
				{
					PL.push_back(pt);
					count++;
				}
				pt.charge=0;
				pt.scan=-1;
			}
		}
	}
	fclose(fp);	
	return count;
}

void PepXML::GetCutoff()
{
	size_t i,sg=FDRTable.size();
	for(i=0;i<sg;i++)
	{
		if(FDRTable[i]>=0.01-1e-6) break;
	}
	if(i<sg)Min_Prob=ErrTable[i];
	else Min_Prob=0;
}