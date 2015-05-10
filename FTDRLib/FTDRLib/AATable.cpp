// AATable.cpp: implementation of the AATable class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "AATable.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
/*charge factor and present frequency
0.069353	0.0200603911755608   
0.078317	0.0369940135537482	
0.068617 	0.0294690644897735
0.048258 	0.0460821392508482
0.067154	0.0426628212053486
0.054213 	0.0288580291114536
0.015276	0.0472636409124965
0.096496	0.0619419465453503
0.059226 	0.0612767540157279
0.041846 	0.0290464650566378	
0.053269	0.0210617795374163
0.039546 	0.0532490794139565
0.059349 	0.26403409842431
0.066427 	0.0351093089604218
0.023801 	0.0605363062979564
0.022930	0.223693977232622
0.040044	0.0672565796079053
0.053549	0.328853420136695
0.030707  	0.0459402814573535
0.011522 	0.0648487956049268
*/
AATable::AATable()
{

	int c='G'-'A';
    table[c].massm=57.02146; table[c].massa=57.0519;
	table[c].rc=-0.9f;table[c].rcnt=5.0f;//	    0.069353  
	table[c].ChC=0.0200603911755608;  
	//---------------------------G---C2H5NO2---------------//
	
	//       C                  H                 N
	table[c].EleNum[0]=2;table[c].EleNum[1]=3;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;
	
	c='A'-'A';
	table[c].massm=71.03711; table[c].massa=71.0788;	 
	table[c].rc=0.8;table[c].rcnt=-1.5;// 	0.078317 
	table[c].ChC=0.0369940135537482;   
	///-----------A----C3H7NO2-------------------------//
	//       C                   H                 N
	table[c].EleNum[0]=3;table[c].EleNum[1]=5;table[c].EleNum[2]=1;
	//         O                  S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;


	c='S'-'A';
	table[c].massm=87.03203; table[c].massa=87.0782;
	table[c].rc=-0.8f;table[c].rcnt=5.0;// 	0.068617 
	table[c].ChC=0.0294690644897735;  
	//---------------------------S---C3H7NO3--------------//
	//       C                  H                 N
	table[c].EleNum[0]=3;table[c].EleNum[1]=5;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=2;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='P'-'A';
	table[c].massm=97.05276; table[c].massa=97.1167;
	table[c].rc=0.2;table[c].rcnt=4.0;// 	    0.048258 
	table[c].ChC=0.0460821392508482;
	//---------------------------P---C5H9NO2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=5;table[c].EleNum[1]=7;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='V'-'A';
	table[c].massm=99.06841; table[c].massa=99.1326;
	table[c].rc=5.0;table[c].rcnt=-5.5;//  	0.067154 
	table[c].ChC=0.0426628212053486;
	//---------------------------V---C5H11NO2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=5;table[c].EleNum[1]=9;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='T'-'A';
	table[c].massm=101.04768;table[c].massa=101.1051; 
	table[c].rc=0.4;table[c].rcnt=5.0;//  	0.054213 
	table[c].ChC=0.0288580291114536;
	//---------------------------T---C4H9NO3---------------//
	//       C                  H                 N
	table[c].EleNum[0]=4;table[c].EleNum[1]=7;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=2;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='C'-'A';
	table[c].massm=103.00919;table[c].massa=103.1388;
	table[c].rc=-0.8;table[c].rcnt=4.0;// 	0.015276 
	table[c].ChC=0.0472636409124965;
	//---------------------------C---C3H7NO2S----------------//
	//       C                   H                 N
	table[c].EleNum[0]=3;table[c].EleNum[1]=5;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=1;table[c].EleNum[5]=0;	

	c='L'-'A';
	table[c].massm=113.08406;table[c].massa=113.1594;
	table[c].rc=9.6f;table[c].rcnt=-9.0;// 	0.096496  
	table[c].ChC=0.0619419465453503;
	//---------------------------L---C6H13NO2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=6;table[c].EleNum[1]=11;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='I'-'A';
	table[c].massm=113.08406;table[c].massa=113.1594;
	table[c].rc=8.4;table[c].rcnt=-8.0;// 	0.059226 
	table[c].ChC=0.0612767540157279;
	//---------------------------I---C6H13NO2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=6;table[c].EleNum[1]=11;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;


	c='X'-'A';
	table[c].massm=113.08406;table[c].massa=113.1594;
	table[c].rc=2.2075;table[c].rcnt=0.8228;//	0.000098  //wight avg all
	table[c].ChC=0.0476769039815862;
	//---------------------------X---L or I---------------//
	//       C                  H                 N
	table[c].EleNum[0]=6;table[c].EleNum[1]=11;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='N'-'A';
	table[c].massm=114.04293;table[c].massa=114.1038;
	table[c].rc=-1.2;table[c].rcnt=5.0;//	0.041846  
	table[c].ChC=0.0290464650566378;
	//---------------------------N---C4H8N2O3---------------//
	//       C                  H                 N
	table[c].EleNum[0]=4;table[c].EleNum[1]=6;table[c].EleNum[2]=2;
	//         O                    S                  P
	table[c].EleNum[3]=2;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='O'-'A';
	table[c].massm=114.07931;table[c].massa=114.1472; 
	table[c].rc=-1.9;table[c].rcnt=4.6;//	0  //use as K  
	table[c].ChC=0.26403409842431;
	//---------------------------O---K--------------//	
	//       C                  H                 N
	table[c].EleNum[0]=4;table[c].EleNum[1]=6;table[c].EleNum[2]=2;
	//         O                    S                  P
	table[c].EleNum[3]=2;table[c].EleNum[4]=0;table[c].EleNum[5]=0;
	
	c='B'-'A';
	table[c].massm=114.53494;table[c].massa=114.5962; 
	table[c].rc=-0.808;table[c].rcnt=7.2402;//	0.000005 //Asn or Asp(D,orN) 
	table[c].ChC=0.0245746550064521f;
	//---------------------------B---D or N---------------//
	//       C                  H                 N
	table[c].EleNum[0]=4;table[c].EleNum[1]=6;table[c].EleNum[2]=2;
	//         O                    S                  P
	table[c].EleNum[3]=2;table[c].EleNum[4]=0;table[c].EleNum[5]=0;


	c='D'-'A';
	table[c].massm=115.02694;table[c].massa=115.0886;
	table[c].rc=-0.5;table[c].rcnt=9.0;//	0.053269 
	table[c].ChC=0.0210617795374163;
	//---------------------------D---C4H7NO4---------------//
	//       C                       H                 N
	table[c].EleNum[0]=4;table[c].EleNum[1]=5;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=3;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='Q'-'A';
	table[c].massm=128.05858;table[c].massa=128.1307;
	table[c].rc=-0.9;table[c].rcnt=1.0;// 	0.039546 
    table[c].ChC=0.0532490794139565;
	//---------------------------Q---C5H10N2O3---------------//
	//       C                  H                 N
	table[c].EleNum[0]=5;table[c].EleNum[1]=8;table[c].EleNum[2]=2;
	//         O                    S                  P
	table[c].EleNum[3]=2;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='K'-'A';
	table[c].massm=128.09496;table[c].massa=128.1741; 
	table[c].rc=-1.9;table[c].rcnt=4.6;// 	0.059349  
	table[c].ChC=0.26403409842431;
	//---------------------------K---C6H14N2O2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=6;table[c].EleNum[1]=12;table[c].EleNum[2]=2;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='Z'-'A';
	table[c].massm=128.55059;table[c].massa=128.6231;	
	table[c].rc=-0.3359;table[c].rcnt=4.761;//	0.000004 ////Q or E 
	table[c].ChC=0.0418785366161028;
	//---------------------------Z--- E or Q--------------//
	//       C                  H                 N
	table[c].EleNum[0]=6;table[c].EleNum[1]=12;table[c].EleNum[2]=2;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='E'-'A';
	table[c].massm=129.04259;table[c].massa=129.1155; 
	table[c].rc=0.0;table[c].rcnt=7.0;// 	0.066427 
	table[c].ChC=0.0351093089604218;
	//---------------------------E---C5H9NO4---------------//
	//       C                  H                 N
	table[c].EleNum[0]=5;table[c].EleNum[1]=7;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=3;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='M'-'A';
	table[c].massm=131.04049;table[c].massa=131.1926; 
	table[c].rc=5.8;table[c].rcnt=-5.5;// 	0.023801  
	table[c].ChC=0.0605363062979564;
	//---------------------------M---C5H11NO2S---------------//
	//       C                  H                 N
	table[c].EleNum[0]=5;table[c].EleNum[1]=9;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=1;table[c].EleNum[5]=0;

	c='H'-'A';
	table[c].massm=137.05891;table[c].massa=137.1411; 
	table[c].rc=-1.3;table[c].rcnt=4.0;// 	0.02293 
	table[c].ChC=0.223693977232622;
	//---------------------------H---C6H9N3O2---------------//
	
	//       C                  H                 N
	table[c].EleNum[0]=6;table[c].EleNum[1]=7;table[c].EleNum[2]=3;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='F'-'A';
	table[c].massm=147.06841;table[c].massa=147.1766;	
	table[c].rc=10.5;table[c].rcnt=-7.0;// 0.040044
	table[c].ChC=0.0672565796079053;   
	//---------------------------F---C9H11NO2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=9;table[c].EleNum[1]=9;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='R'-'A';
	table[c].massm=156.10111;table[c].massa=156.1875;	
	table[c].rc=-1.3;table[c].rcnt=8.0;// 	0.053549 
	table[c].ChC=0.328853420136695;
	//---------------------------R---C6H14N4O2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=6;table[c].EleNum[1]=12;table[c].EleNum[2]=4;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='Y'-'A';
	table[c].massm=163.06333;table[c].massa=163.1760;
	table[c].rc=4.0;table[c].rcnt=-3.0;//	0.030707   
	table[c].ChC=0.0459402814573535;
	//---------------------------Y---C9H11NO3---------------//
	//       C                  H                 N
	table[c].EleNum[0]=9;table[c].EleNum[1]=9;table[c].EleNum[2]=1;
	//         O                    S                  P
	table[c].EleNum[3]=2;table[c].EleNum[4]=0;table[c].EleNum[5]=0;

	c='W'-'A';
	table[c].massm=186.07931;table[c].massa=186.2132;
	table[c].rc=11.0;table[c].rcnt=-4.0;// 0.011522  
	table[c].ChC=0.0648487956049268;
	//---------------------------W---C11H12N2O2---------------//
	//       C                  H                 N
	table[c].EleNum[0]=11;table[c].EleNum[1]=10;table[c].EleNum[2]=2;
	//         O                    S                  P
	table[c].EleNum[3]=1;table[c].EleNum[4]=0;table[c].EleNum[5]=0;
	// the element compose
	
}

AATable::~AATable()
{

}

void AATable::Add2Mass(char AAName, double deltmass,int masstype)
{
	int c=AAName-'A';
	if(c>=0&&c<26)
	{
		if(masstype==MASS_AVG)	table[c].massa+=deltmass;
		else if(masstype==MASS_MON) table[c].massm+=deltmass;		
	}
}

void AATable::SetMass(char AAName,double mass,int masstype)
{
	int c=AAName-'A';
	if(c>=0&&c<26)
	{
		if(masstype==MASS_AVG)	table[c].massa=mass;
		else if(masstype==MASS_MON) table[c].massm=mass;		
	}
}

double AATable::GetMass(char AAName,int masstype)
{
	double mass;
	if(AAName=='n') 
	{
		if(masstype==MASS_AVG)mass=1.00794;
		else if(masstype==MASS_MON) mass=1.007825035;
	}
	else if(AAName=='c') 
	{
		if(masstype==MASS_AVG)return 17.00734;
		else if(masstype==MASS_MON) mass=17.002739665;
	}
	int c=AAName-'A';	
	if(c>=0&&c<26) 
	{
		if(masstype==MASS_AVG)mass=table[c].massa;
		else if(masstype==MASS_MON) mass=table[c].massm;
	}
    return mass;
}

double AATable::GetRC(char AAName)
{
	int c=AAName-'A';
	if(c>=0&&c<26) 	return table[c].rc;
	return 0;
}

double AATable::GetRCNT(char AAName)
{
	int c=AAName-'A';
	if(c>=0&&c<26) 	return table[c].rcnt;
	return 0;
}

double AATable::GetChC(char AAName)
{
	int c=AAName-'A';
	if(c>=0&&c<26) 	return table[c].ChC;
	return 0;
}	

bool AATable::GetEleCom(int *EleComp,char AAName)
{
	int c=AAName-'A';
	if(c>=0&&c<26)
	{
		 for(int i=0;i<6;i++) EleComp[i]=table[c].EleNum[i];
		 return true;
	}
	return false;
}

void AATable::GetEleCom(int *EleComp,CString seq)
{
	int i;
	for(i=0;i<6;i++) EleComp[i]=0;
	int slen=seq.GetLength();
	for(i=0;i<slen;i++)
	{
		char c=seq.GetAt(i)-'A';
		if(c>=0&&c<26)
		{
			for(int k=0;k<6;k++) EleComp[k]+=table[c].EleNum[k];
		}
		//take into account the fixed modification
		if(c=='C'-'A')
		{	
	 
		   EleComp[0]+=3;
		   EleComp[1]+=2;
		   EleComp[2]+=1;
		   EleComp[3]+=1;
		}
	}
	//add a H2O
	EleComp[1]+=2;
	EleComp[3]+=1;
}

AA::AA()
{
	massm=0;
	massa=0;
	rc=0;
	rcnt=0;
	ChC=0;
	for(int i=0;i<6;i++)
	{
		EleNum[i]=0;
	}
}

AA::~AA()
{
}

AA &AA::operator=(const AA &r)
{
	massm=r.massm;
	massa=r.massa;
	rc=r.rc;
	rcnt=r.rcnt;
	ChC=r.ChC;
	for(int i=0;i<6;i++)
	{
		EleNum[i]=r.EleNum[i];
	}
	return *this;
}

AA::AA(const AA &r)
{
	massm=r.massm;
	massa=r.massa;
	rc=r.rc;
	rcnt=r.rcnt;
	ChC=r.ChC;
	for(int i=0;i<6;i++)
	{
		EleNum[i]=r.EleNum[i];
	}
}