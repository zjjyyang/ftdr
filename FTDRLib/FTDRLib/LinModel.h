#pragma once
#include<vector>
using std::vector;
#include "commonddef.h"
#include "GaussFit.h"
class LinModel
{
public:
	double *COFF;
	int PARNUM;
public:	
	double stats[3];
	LinModel(void);
	~LinModel(void);
	LinModel(const LinModel &cp);
	LinModel &operator=(const LinModel &cp);
	int scanSize(vector<CalData> &DataList);
	bool ModelTrain(vector<CalData> &DataList,vector<size_t> &ParSel);
	double Predict(CalData &dt,vector<size_t> &ParSel);
	void Predict(CalData &dt,vector<double> &ppmD,vector<size_t> &ParSel);	
	void Predict(CalData &dt,vector<double> &ppmD,vector<double> &weight,vector<size_t> &ParSel);
	double Predict(double *fea,int FEATURE_NUM);
};
