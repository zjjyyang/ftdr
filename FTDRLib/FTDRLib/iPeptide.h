// iPeptide.h: interface for the iPeptide class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#define COLOR_NON 0x00
#define COLOR_RED 0x01
#define COLOR_BLACK 0x02

class iPeptide  
{
public:
	 double Observed;
	 double CalcMass;
	 double ObsrvMass;
	 double Expect;
	 double dm;
	 int miss;
	 int score;
	 int rank;	
	 char color;
	 bool IsBold;
	 bool IsChecked;
	 int ScanNumber;
	 CString fName;
	 CString Sequence;
	 CString Modif; 
public:
	 iPeptide();
	 ~iPeptide();	 
	 iPeptide &operator=(iPeptide &r);
	 void Clear();	
};

