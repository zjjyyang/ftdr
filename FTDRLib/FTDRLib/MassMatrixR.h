#pragma once
#include "commonddef.h"
#include <vector>
using std::vector;
class MassMatrixR
{
public:
	MassMatrixR(void);
	~MassMatrixR(void);
	int LoadData(char *fname,vector<Pre_calData> &PL);
	char *seekTable(int n,char *buf);
	bool freadline(FILE *fp,char *buf,int maxcount);
};
