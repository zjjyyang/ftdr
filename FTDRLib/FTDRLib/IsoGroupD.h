#pragma once
#include "mPeak.h"
#include "stdio.h"
#include <gsl/gsl_vector.h>
#include <vector>
#include <string>
using std::vector;
using std::string;

//for IsoGroupD
///////////////////////
#define NO_PEAKS 0
#define NO_GOOD_CLUSTER 1
#define CPEAK_FAILUR 2
#define FIND_ONE 3
////////////////////
#define PRE_DEF_WT 0.2
#define LOW_FACTOR 0.5
#define GD_HIGH 0.8
#define GD_LOW_CT 0.3
#define FACTOR_MUL 1.5
#define MIN_PEAK_INT 1
#define NAN 1e30
#define ABS_SIG_CUT 0.01
#define F_K 1000
#define SIGNAL_FACTOR 1e6
///////////////////for the parent ion dectection
class IsoLocal
{
public:
	int UnitIdx;
	int UnitLocal;
	IsoLocal();
	~IsoLocal();
	IsoLocal(const IsoLocal &cp);
	IsoLocal& operator =(const IsoLocal &cp);

};

class sPeak
{
public:
	double mz;
	double dInt;
	double base;
	double noise;
	bool IsUsed;
	int Global_IDX;
	vector<IsoLocal> LocalList;
	int count;
	sPeak();
	~sPeak();
	sPeak(const sPeak &r);
	sPeak(const mPeak &r);
	sPeak &operator=(const sPeak &r);	
	sPeak &operator=(mPeak &r);	
	bool IsInLocalList(int idx);
};

class refPeak
{
public:
	double mz;
	double dInt;
	int clusterid;
	int localid;
	int IsoNum;
	double S_N;
	int rank;
	double weight;
	refPeak();
	~refPeak();
	refPeak(const refPeak &cp);
	refPeak &operator=(const refPeak &cp);
	refPeak &operator=(const sPeak &cp);
};

typedef vector<sPeak> MS_PKL;

class IsoUnit
{
public:
	int charge;
	vector<int> peak_idx;
	int type;
	double weight;
	IsoUnit();
	~IsoUnit();
	IsoUnit(const IsoUnit& cp);
	IsoUnit& operator=(const IsoUnit& cp);
	bool CheckFactor(IsoUnit &cs);
	bool CheckCover(IsoUnit &cs);
	void AdjustCH(double mass);
};

class FitVector
{
public:
	double Exp_Int;
	vector<size_t> IDX;      //]must be the same size
	vector<double> IsoTL;// ]
	FitVector();
	~FitVector();
	FitVector(const FitVector &r);
	FitVector &operator=(const FitVector &r);
	int FindIdx(int idx);
};

typedef vector<FitVector> fParList;

class IsoGroup
{
public:
	MS_PKL gPeaks;
	vector<IsoUnit> UList;	
	IsoGroup();
	~IsoGroup();
	IsoGroup(const IsoGroup& cp);
	IsoGroup& operator=(const IsoGroup& cp);
	void Initial(sPeak &st);
	bool AddPeaks(sPeak &st,double ppmDM);
	void OutPut();
	void SortUList();
	//void PostSplit();
	void UpdateWeight(double Pre_def_w);
	double GetUnitGD(int ui);
	double GetCOS(double *EInt,double *TInt,int num);
	bool IsANew(double mass, double Iso1,double Iso2,int idx);
	//void GetIsoDis(double mass,double IsoDisT[MAX_ISO]);
	//void GetIsoDis(double mass,double *IsoDisT,int max_iso);
	//double GetIsoDis(double mass,int iso_no);
	void SplitCH();	
	bool InitialPars(vector<refPeak> &ACIs,fParList &PL);
	double Fitting(double gd_cut,int Min_Iso_Num,vector<CalData> &XICPL);
	void ReInitialByLS(vector<refPeak> &ACIs,fParList &PL);
	void GetCharge(vector<int> &ch);

};

double Res_f(const gsl_vector *v, void *params);
void Res_df(const gsl_vector *v, void *params,gsl_vector *df);
void Res_fdf(const gsl_vector *v, void *params,double *f,gsl_vector *df);
void Res_feps(const gsl_vector *v, void *params,double *eps);
