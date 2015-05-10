#pragma once
#include "svm.h"
#include<vector>
#include "commonddef.h"
using std::vector;
#define TRY_PAR_NUM 4
const static int CTry[TRY_PAR_NUM]={64,128,256,512};
class SVMModel
{
private:
	struct svm_model *SVM_MODEL;
	vector<double> feature_max;
	vector<double> feature_min;
	//must be keep for the future use of the model
	struct svm_node *svm_x_space;
	struct svm_parameter svm_param;		// set by parse_command_line
	struct svm_problem svm_prob;		// set by read_problem
	//double C_par;
private:
	int scanSize(vector<CalData> &DataList);	
	double CalRsquare(int nr_fold);
public:
	double Res_mean;
	double Res_std;
	double gd;
	SVMModel(void);
	~SVMModel(void);	
	void Scale(int Item_Num,int FEATURE_NUM,struct svm_node *svm_x_space);	
	void Scale_P(svm_node *x);
	void Scale_P(svm_node *x,int FEATURE_NUM);
	bool Model_Train(vector<CalData> &DataList,vector<size_t> &ParSel);
	bool RMSD(int n,int FEATURE_NUM,struct svm_node *svm_x_space,svm_problem *svm_prob);
	void Predict(CalData &dt,vector<double> &ppmD,vector<size_t> &ParSel);
	double Predict(CalData &dt,vector<size_t> &ParSel);	
	void Predict(CalData &dt,vector<double> &ppmD,vector<double> &weight,vector<size_t> &ParSel);
	double Predict(double *fea,int FEATURE_NUM);
	bool OutputModel(const char *fname);	
};
