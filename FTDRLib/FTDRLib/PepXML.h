#pragma once

#include "commonddef.h"
#include <vector>
using std::vector;

class PepXML
{
private:
	double Min_Prob;
	void RemTail(char *buf);
	vector<double> ErrTable;
	vector<double> FDRTable;
public:
	PepXML(void);
	~PepXML(void);
	int conf_level;	
	int LoadData(char *fname,vector<Pre_calData> &PL);
	void GetCutoff();
};
