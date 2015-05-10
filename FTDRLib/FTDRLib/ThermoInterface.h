// -*- mode: c++ -*-


/*
    File: ThermoInterface.h
    Description: Encapsulation for Thermo Xcalibur interface.
    Date: July 25, 2007

    Copyright (C) 2007 Natalie Tasman, ISB Seattle


    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/


#pragma once

#include <string>
#include "xrawfile2_2_0.h"
#include "InstrumentInterface.h"
#include "FilterLine.h"
#include "SVMModel.h"
#include "LinModel.h"
#include "commonddef.h"


typedef struct _datapeak
{
	double dMass;
	double dIntensity;
} DataPeak;


class ThermoInterface : public InstrumentInterface {
private:
	// COleDispatchDriver object instance to xrawfile2 dll
	
	// this is the smart pointer we'll use as the entry point;
	// it may actually be initialized as IXRawfile2 or IXRawfile3
	// depending on runtime-detected installed Xcalibur DLL version
	XRAWFILE2Lib::IXRawfile2Ptr xrawfile2_;
	int IXRawfileVersion_; // which IXRawfile interface was initialized?

	int msControllerType_;
	long firstScanNumber_;
	long lastScanNumber_;
	bool firstTime_;
	bool getPrecursorInfo(Scan& scan, long scanNumber, FilterLine& filterLine);
	bool forcePrecursorFromFilter_;
	bool getStatusData(Scan *cs);
	bool getTailerData(Scan *cs);

public:
	int getPreInfoCount_;
	int filterLineCount_;
	int oldAPICount_;

public:
	long curScanNum_;

public:
	ThermoInterface(void);
	~ThermoInterface(void);

	virtual bool initInterface(void);
	virtual bool setInputFile(const std::string& fileName);
	virtual void setCentroiding(bool centroid);
	virtual void setDeisotoping(bool deisotope);
	virtual void setCompression(bool compression);
	virtual void forcePrecursorFromFilter(bool mode);
	virtual void setVerbose(bool verbose);
	virtual void setShotgunFragmentation(bool sf) {}
	virtual void setLockspray(bool ls) {}
	virtual Scan* getScan(void);
	//added from 2011.3.20, by zhangjy
	//get a scan data for the given scannumber
	Scan* getScan(long scannumber);
	long getFirstScanNumber();
	long GetLastScanNumber();
	int MSLevel(long ScanNumber);
	int MSLevel();
	bool SkipScan();
	//for ms 1 calibrate 
public:
	SVMModel *CalModel;
	bool IsCalibrate;
	void SortMerge(vector<double> &tmz,vector<double> &tint);
	double GetSampleDm(double mass);
	void SeekToScan(long scan);
	LinModel LM;
	int MType;
	double GetCalibrateMZ(double *fea,int FEA_NUM);
public:
	vector<size_t> statusIDX;
	size_t trailerIDX[2];	
	void Convert(vector<size_t> &ColParIdx,vector<size_t> &SelParIdx,size_t tIDX[2]);
	void CalibeateMS_SVM(Scan *curScan);
	void CalibeateMS_LIN(Scan *curScan);
public:
	double MzIntModel[4];
	//intrval=MzIntModel[2]*mz^2+MzIntModel[1]*mz+MzIntModel[0];
	//R2=MzIntModel[3]
	bool PreSuryMZINT();
public:
	int TFindH(DataPeak *PKL,int b,int e,double fmz);
	int TFindL(DataPeak *PKL,int b,int e,double fmz);	
	int LocatePmz(DataPeak *PKL,int n,double begin,double end,double pmz);	
	int ReselectPIons(DataPeak *PKL,int dpnum,double curpmz,int ch);
public:
	vector<double> CalPmz;
	vector<long> scanmap;
	long totalpmzscan;
	long seekScan(long scannum,size_t a, size_t b);
	//bool IsUseCpmz;
public:
	void ReInitial();
public:
	int countMS2();
	int countMS1();
	void countMS(int CT[2]);
};


