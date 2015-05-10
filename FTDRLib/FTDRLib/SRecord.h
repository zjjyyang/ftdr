// SRecord.h: interface for the SRecord class.
//

#pragma once
class SRecord  
{
public:
	SRecord();
	~SRecord();
	int RankSp;
	float MH;
	float detCn;
	float XCorr;
	float Sp;
	int mIons;
	int tIons;
	CString Ref;
	CString Seq;
	SRecord &operator =(SRecord &r);
};

