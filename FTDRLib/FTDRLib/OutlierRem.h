#pragma once
#include <vector>
using std::vector;
#define INF 1e-16

class OutlierRem
{
private:
	double _K;
	vector<int> _w;
	vector<double>_x;
	double _mean();
	double _std(double m);
	bool _updatew();
	void initialw(vector<double> *data);
public:	
	double m_new;
	double std_new;
	OutlierRem(void);
	~OutlierRem(void);
	bool oRemove(vector<double> &x);
	bool oRemoveGTest(vector<double> &x);
	bool oRemoveGTestW(vector<double> &x,vector<double> &w);
	bool oRemoveChauCt(vector<double> &x);
	bool ORemMixModel(vector<double> &x);
	bool oRemoveGTest(double *y, int N);
	bool oRemoveChauCt(double *y, int N);
	int count(double *x,int n,double min_v, double max_v);
	double gsl_sum(double *x,int n);
	bool ORemMixModel(double *x,int n);	
	double GetMean();
	double GetStd();	
	double Mean(vector<double> &x);
	double Std(double m,vector<double> &x);
	double MeanW(vector<double> &x,vector<double> &w);
	double StdW(double m,vector<double> &x,vector<double> &w);
};
