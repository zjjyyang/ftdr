#pragma once
#include <gsl/gsl_histogram.h>
#include <vector>
using std::vector;

#define C0 0.39894228040143

#define PRE_MC 50
#define ADD_CT 20

class myHistogram
{
public:
	double *range;
	double *bin;
	int n;
	double mean;
	double std;
	CString rawname;
	myHistogram();
	myHistogram(const myHistogram &r);
	myHistogram(gsl_histogram *r);
	~myHistogram();
	myHistogram &operator=(const myHistogram &r);
	bool Convert(double *res,int n);
	void Convert(gsl_histogram *r,int num);
	bool Convert(vector<double> &x);
	double GetMaxY();
	double GetMinX();
	double GetMaxX();
	double GetProbGS(double x);
	bool alloc(int num);
};

class HistList
{
private:
	myHistogram *HPL;
	int Mcount;
	int Acount;
public:
	HistList(void);
	HistList(int n);
	HistList(HistList &r);
	~HistList(void);
	int size();
	void reshape();
	void clear();
	bool at(int pos,myHistogram &mt);
    bool AddHist(gsl_histogram *r,double m,double sd,CString raw);
	bool AddHist(double *res,int n,double m,double sd,CString raw);
	void push_back(const myHistogram mt);
	HistList &operator=(HistList &r);	
};
