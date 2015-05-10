// mPeak.h: interface for the mPeak class.

#pragma once
#include "commonddef.h"

class mPeak  
{
public:
	mPeak();
	~mPeak();
	double dMass;
	double dIntensity;
	double rInt;	
	double baseLine;
	double noise;
	//end
	mPeak(const mPeak &r);
	mPeak &operator=(const mPeak &r);	
	//this function must be safe call by others
	void Write2File(FILE *fp);//the size is sizeof(double)*5, about 40 generally
	//fp must be read and write with bineary model
	void ReadFromFile(FILE *fp);
};

