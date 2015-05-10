#pragma once
#define MASS_AVG 0
#define MASS_MON 1
///////////////氨基酸残基的定义/////////////////////////
///ref:
//Title:An Improved Model for Prediction of RetentionTimes of Tryptic Peptides in Ion Pair Reversed-phase HPLC
//Joural:Molecular & Cellular Proteomics 3.9 page:912 ,2004
//Author:O.V.Krokhin R.Craig V.Spicer W.Ens, K.G.Standing, R.C.Beavisand J.A.Wilkins
// the AA percent come frome sp database release version 48.7

class AA
{
public:
	AA();
	~AA();
	double massm;
	double massa;
	double rc;
	double rcnt;
	double ChC;
	int EleNum[6];
	AA(const AA &r);
	AA &operator=(const AA &r);
};

class AATable  
{
private:
	AA table[26];
public:
	AATable();
	~AATable();
	void Add2Mass(char AAName, double deltmass,int masstype);
	void SetMass(char AAName,double mass,int masstype);
	double GetMass(char AAName,int masstype);
	double GetRC(char AAName);
	double GetRCNT(char AAName);
	double GetChC(char AAName);
	bool GetEleCom(int *EleComp,char AAName);
	void GetEleCom(int *EleComp,CString seq);
};



