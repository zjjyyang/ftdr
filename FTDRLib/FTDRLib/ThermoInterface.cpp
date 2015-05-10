// -*- mode: c++ -*-


/*
File: ThermoInterface.cpp
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


#include <iostream>
#include <iomanip>

#include "stdafx.h"
#include "ThermoInterface.h"
#include "Scan.h"
#include "MSUtilities.h"
#include "math.h"
#include "CalScan.h"

using namespace std;
#include "OTrace.h"

#include "GaussFit.h"

ThermoInterface::ThermoInterface(void)
{

	// InstrumentInterface members
	totalNumScans_ = -1;
	curScanNum_ = -1;
	firstScanNumber_ = -1;
	lastScanNumber_ = -1;
	doCompression_ =false;
	doCentroid_=false;
	doDeisotope_=false;
	forcePrecursorFromFilter_=false;
	verbose_=false;
	accurateMasses_=0;
	inaccurateMasses_=0;
	// store counts for charges up to +15
	chargeCounts_.clear();
	chargeCounts_.resize(1, 0);
	instrumentInfo_.manufacturer_ = THERMO;
	instrumentInfo_.acquisitionSoftware_ = XCALIBUR;

	// ThermoInterface members
	xrawfile2_ = NULL; // smartptr will be initialized in initInterface
	IXRawfileVersion_ = -1;
	msControllerType_ = 0;
	startTimeInSec_ = -1;
	endTimeInSec_ = -1;

	firstTime_ = true;

	getPreInfoCount_ = 0;
	filterLineCount_ = 0;
	oldAPICount_ = 0;
	///// addtional by zhangjy 2011.9.5
	CalModel=NULL;
	IsCalibrate=false;

	//a default model is for Normal Orbitrap 
	MzIntModel[0]=1.12e-009;
	MzIntModel[1]=2.204e-006;
	MzIntModel[2]=-0.0003366;
	MzIntModel[3]=1.0;

	MType=MODEL_SVM;

	trailerIDX[0]=0;
	trailerIDX[1]=0;

	totalpmzscan=0;

	//IsUseCpmz=false;
}

ThermoInterface::~ThermoInterface(void)
{
	xrawfile2_->Close();
}

void ThermoInterface::ReInitial()
{
	curScanNum_ = firstScanNumber_;
	firstTime_ = true;
}

bool ThermoInterface::initInterface(void) 
{

	// temp smartptrs; xrawfile2_ will be set to highest interface that succeeds
	XRAWFILE2Lib::IXRawfile2Ptr xrawfile2IXRawfile2_ = NULL;
	XRAWFILE2Lib::IXRawfile3Ptr xrawfile2IXRawfile3_ = NULL;

	// first, init MFC or return false
	//if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	//{
	//	cerr << _T("Fatal Error: MFC initialization failed") << endl;
	//	return false;
	//}

		// Initializes the COM library on the current thread 
	// and identifies the concurrency model as single-thread 
	// apartment (STA)
	CoInitialize( NULL );		

	// first try to use the IXRawfile3 interface, introduced in Xcalibur 2.0
	HRESULT hr = xrawfile2IXRawfile3_.CreateInstance("XRawfile.XRawfile.1");
	if (!FAILED(hr)) {
		// we can use the IXRawfile3 interface
		xrawfile2_ = xrawfile2IXRawfile3_;
		IXRawfileVersion_  = 3;
		//cout << "Xcalibur 2.0 interface initialized." << endl;
		printf("Xcalibur 2.0 interface initialized.\n");
	}
	else {
		//cerr << "Unable to initialize Xcalibur 2.0 interface; falling back to 1.4 interface" << endl;
		printf("Unable to initialize Xcalibur 2.0 interface; falling back to 1.4 interface\n");
		// next, try to use the IXRawfile2 interface, found in Xcalibur 1.4
		HRESULT hr = xrawfile2IXRawfile2_.CreateInstance("XRawfile.XRawfile.1");
		if (!FAILED(hr)) {
			// we can use the IXRawfile2 interface
			xrawfile2_ = xrawfile2IXRawfile2_;
			IXRawfileVersion_ = 2;
			//cout << "Xcalibur 1.4 interface initialized." << endl;
			printf("Xcalibur 1.4 interface initialized.\n");
		}
		else {
			//cerr << "Unable to initialize XCalibur 1.4 interface; check installation" << endl;
			//cerr << "hint: try running the command \"regsvr32 C:\\<path_to_Xcalibur_dll>\\XRawfile2.dll\"" << endl;
			printf("Unable to initialize XCalibur 1.4 interface; check installation\n");
			printf("hint: try running the command \"regsvr32 C:\\<path_to_Xcalibur_dll>\\XRawfile2.dll\\n");
			return false;
		}
	}

	// otherwise, ok!
	return true;
}

bool ThermoInterface::setInputFile(const string& filename) {

	// try to open raw file or die
	HRESULT hr = xrawfile2_->Open(filename.c_str());
	if (hr != ERROR_SUCCESS) {

		if (hr == ERROR_PATH_NOT_FOUND) {
			//cerr << filename << ": path not found" << endl;
			printf("%s : path not found\n",filename.c_str());
		}
		if (hr == ERROR_FILE_NOT_FOUND) {
			//cerr << filename << ": file not found" << endl;
			printf("%s : file not found\n",filename.c_str());
		}

		//cerr << "Cannot open file " << filename << endl;
		printf("Cannot open file %s\n",filename.c_str());
		//cerr << "(If you get this error with a valid filename," << endl;
		printf("(If you get this error with a valid filename,\n");
		//cerr << " you may be using an older version of Xcalibur libraries" << endl;	
		printf(" you may be using an older version of Xcalibur libraries\n");
		//cerr << " with a file created by a newer version of Xcalibur.)" << endl;	
		printf(" with a file created by a newer version of Xcalibur.)\n");
		
		return false;
	}


	//cout << "(Thermo lib opened file " << filename << ")" << endl;
	printf("(Thermo lib opened file %s\n",filename.c_str());

	// test the file format version number
	long fileVersionNumber = -1;
	xrawfile2_->GetVersionNumber(&fileVersionNumber);
	//cout << "file version is " << fileVersionNumber << ", interface version is ";
	printf("file version is %d, interface version is ",fileVersionNumber);
	if (IXRawfileVersion_ == 2) {
			//cout << "1.4" << endl;
		printf("1.4\n");
	}
	else if (IXRawfileVersion_ == 3) {
			//cout << "2.0 or greater" << endl;
		printf("2.0 or greater\n");
	}
	

	xrawfile2_->SetCurrentController(msControllerType_, 1);
	// Get the total number of scans

	xrawfile2_->GetFirstSpectrumNumber(&firstScanNumber_);
	xrawfile2_->GetLastSpectrumNumber(&lastScanNumber_);
	//cout << "file should contain scan numbers " << firstScanNumber_ << " through " << lastScanNumber_ << endl;
	printf("file should contain scan numbers %d through %d\n",firstScanNumber_,lastScanNumber_);

	curScanNum_ = firstScanNumber_;

	totalNumScans_ = (lastScanNumber_ - firstScanNumber_) + 1;

	// get the start and the end time of the run
	double startTime, endTime;
	xrawfile2_->GetStartTime(&startTime);
	xrawfile2_->GetEndTime(&endTime);
	// convert from minutes to seconds
	startTimeInSec_ = 60.0 * startTime;
	endTimeInSec_ = 60.0 * endTime;

	// get the instrument model
	BSTR bstrInstModel=NULL;
	xrawfile2_->GetInstModel(&bstrInstModel);
	string instModel = toUpper(convertBstrToString(bstrInstModel));
	SysFreeString(bstrInstModel);
	if (instModel == "LTQ") {
		instrumentInfo_.instrumentModel_ = LTQ;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "LTQ XL") {
		instrumentInfo_.instrumentModel_ = LTQ_XL;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "LTQ FT") {
		instrumentInfo_.instrumentModel_ = LTQ_FT;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "LTQ FT ULTRA") {
		instrumentInfo_.instrumentModel_ = LTQ_FT_ULTRA;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "LTQ ORBITRAP") {
		instrumentInfo_.instrumentModel_ = LTQ_ORBITRAP;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "LTQ ORBITRAP DISCOVERY") {
		instrumentInfo_.instrumentModel_ = LTQ_ORBITRAP_DISCOVERY;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "LTQ ORBITRAP XL") {
		instrumentInfo_.instrumentModel_ = LTQ_ORBITRAP_XL;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "LXQ") {
		instrumentInfo_.instrumentModel_ = LXQ;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "TSQ QUANTUM ACCESS") {
		instrumentInfo_.instrumentModel_ = TSQ_QUANTUM_ACCESS;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "TSQ QUANTUM") {
		instrumentInfo_.instrumentModel_ = TSQ_QUANTUM;
		instrumentInfo_.manufacturer_ = THERMO_SCIENTIFIC;
	}
	else if (instModel == "LCQ") {
		instrumentInfo_.instrumentModel_ = LCQ;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN; // TODO: verify "finnigan"
	}
	else if (instModel == "LCQ ADVANTAGE") {
		instrumentInfo_.instrumentModel_ = LCQ_ADVANTAGE;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "LCQ CLASSIC") {
		instrumentInfo_.instrumentModel_ = LCQ_CLASSIC;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "LCQ DECA") {
		instrumentInfo_.instrumentModel_ = LCQ_DECA;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}	
	else if (instModel == "LCQ DECA XP") {
		instrumentInfo_.instrumentModel_ = LCQ_DECA_XP;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "LCQ DECA XP PLUS") {
		instrumentInfo_.instrumentModel_ = LCQ_DECA_XP_PLUS;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "LCQ FLEET") {
		instrumentInfo_.instrumentModel_ = LCQ_FLEET;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "TSQ VANTAGE STANDARD") {
		instrumentInfo_.instrumentModel_ = TSQ_VANTAGE_STANDARD;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "LTQ VELOS") {
		instrumentInfo_.instrumentModel_ = LTQ_VELOS;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else if (instModel == "LTQ ORBITRAP VELOS") {
		instrumentInfo_.instrumentModel_ = LTQ_ORBITRAP_VELOS;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}
	else {
		instrumentInfo_.instrumentModel_ = INSTRUMENTMODEL_UNDEF;
		instrumentInfo_.manufacturer_ = THERMO_FINNIGAN;
	}


	// get instrument name
	BSTR bstrInstName=NULL;
	xrawfile2_->GetInstName(&bstrInstName);
	//cout << convertBstrToString(bstrInstName) << endl;
	instrumentInfo_.instrumentName_ = convertBstrToString(bstrInstName);
	SysFreeString(bstrInstName);

	// TODO: is this the Xcalibur or inst version?
	// get acquisition software version
	BSTR bstrAcquSoftwareVersion=NULL;
	xrawfile2_->GetInstSoftwareVersion(&bstrAcquSoftwareVersion);
	instrumentInfo_.acquisitionSoftwareVersion_ = convertBstrToString(bstrAcquSoftwareVersion);
	SysFreeString(bstrAcquSoftwareVersion);

	// get instrument hardware version
	BSTR bstrInstHardwareVersion=NULL;
	xrawfile2_->GetInstHardwareVersion(&bstrInstHardwareVersion);
	//cout << convertBstrToString(bstrInstHardwareVersion) << endl;
	instrumentInfo_.instrumentHardwareVersion_ = convertBstrToString(bstrInstHardwareVersion);
	SysFreeString(bstrInstHardwareVersion);

	// get instrument description
	BSTR bstrInstDescription=NULL;
	xrawfile2_->GetInstrumentDescription(&bstrInstDescription);
	//cout << convertBstrToString(bstrInstDescription) << endl;
	SysFreeString(bstrInstDescription);

	// get instrument serial number
	BSTR bstrInstSerialNumber=NULL;
	xrawfile2_->GetInstSerialNumber(&bstrInstSerialNumber);
	//cout << convertBstrToString(bstrInstSerialNumber) << endl;
	SysFreeString(bstrInstSerialNumber);

	// get instrument ion source, analyzer 
	// assuming first scan's filter line info is true for all scans
	// TODO: This is known to be wrong in LTQ-FT case (both FT and ion trap)
	BSTR bstrFilter = NULL;
	xrawfile2_->GetFilterForScanNum(1, &bstrFilter);
	string filterString = convertBstrToString(bstrFilter);
	SysFreeString(bstrFilter);
	FilterLine filterLine;
	if (!filterLine.parse(filterString)) 
	{
		//cerr << "error parsing filter line. exiting." << endl;
		printf( "error parsing filter line. exiting.\n");
		//cerr << "line: " << filterString << endl;
		printf( "line: %s\n",filterString.c_str());
		exit(-1);
	}

	// TODO: deal with muliple analyzers, like FT/IT
	MSAnalyzerType analyzer = filterLine.analyzer_;
	// TODO: hack! do this properly with library somehow?
	if (analyzer == ANALYZER_UNDEF) {
		// set default
		analyzer = ITMS;
	}
	instrumentInfo_.analyzerList_.push_back(analyzer);

	MSIonizationType ionization = filterLine.ionizationMode_;
	// TODO: hack! determine instrument info better
	if (ionization == IONIZATION_UNDEF) {
		// set default
		ionization = NSI;
	}
	instrumentInfo_.ionSource_ = ionization;

	// get time from file
	DATE date=NULL;
	xrawfile2_->GetCreationDate(&date);
	//instrumentInfo_.runTime_ = ""; // FIX
	//CTime c(date);
	COleDateTime dateTime(date);
	string dateStr=dateTime.Format("%Y-%m-%dT%H:%M:%S");
	timeStamp_ = dateStr;

	return true;
}

void ThermoInterface::setVerbose(bool verbose) {
	verbose_ = verbose;
}

void ThermoInterface::setCentroiding(bool centroid) {
	doCentroid_ = centroid;
}


void ThermoInterface::setDeisotoping(bool deisotope) {
	doDeisotope_ = deisotope;
}


void ThermoInterface::forcePrecursorFromFilter(bool force) {
	forcePrecursorFromFilter_ = force;
}

void ThermoInterface::setCompression(bool compression) {
	doCompression_ = compression;
}

//get a scan data for the given scannumber
Scan* ThermoInterface::getScan(long scannumber)
{
	if (scannumber > lastScanNumber_) return NULL;
	Scan* curScan = new Scan();
	curScan->isThermo_ = true;	
	//curScan->curScanNum=scannumber;
	//// test the "scan event" call
	//// gives a more through scan filter line
	//BSTR bstrScanEvent = NULL;
	//xrawfile2_->GetScanEventForScanNum(curScanNum_, &bstrScanEvent);
	//SysFreeString(bstrScanEvent);

	// Get the scan filter
	// (ex: "ITMS + c NSI Full ms [ 300.00-2000.00]")

	BSTR bstrFilter = NULL;
	xrawfile2_->GetFilterForScanNum(scannumber, &bstrFilter);
	curScan->thermoFilterLine_ = convertBstrToString(bstrFilter);
	SysFreeString(bstrFilter);
	FilterLine filterLine;
	if (!filterLine.parse(curScan->thermoFilterLine_)) 
	{
		//cerr << "error parsing filter line. exiting." << endl;
		printf("error parsing filter line. exiting.\n");
		//cerr << "line: " << curScan->thermoFilterLine_ << endl;
		printf("line: %d\n",curScan->thermoFilterLine_);
		exit(-1);
	}

	// we should now have:
	// msLevel
	// polarity
	// scanType (full, zoom, etc)
	// ionization type
	// analyzer
	// scan range
	// parent mass and CID energy, if MS >=2 

	// record msLevel from filter line
	curScan->msLevel_ = filterLine.msLevel_;

	curScan->curScanNum=scannumber;//added by zhangjy on 2011.1.4

	// record polarity from filter line
	curScan->polarity_ = filterLine.polarity_;

	// record analyzer from filter line
	curScan->analyzer_ = filterLine.analyzer_;

	// record ionization from filter line
	curScan->ionization_ = filterLine.ionizationMode_;

	// record scan type from filter line
	// (zoom, full, srm, etc)
	curScan->scanType_ = filterLine.scanType_;

	// record activation (CID, etc)
	// check FilterLine: this may be default set to CID now
	if (filterLine.activationMethod_ != ACTIVATION_UNDEF) {
		curScan->activation_ = filterLine.activationMethod_;
	}
	else {
		curScan->activation_ = CID;
	}

	// record scan ranges from filter line
	// Note: SRM fills this with range of q3 transtion lists
	if (curScan->scanType_ != SRM) {
		curScan->startMZ_ = filterLine.scanRangeMin_[ filterLine.scanRangeMin_.size() - 1 ];
		curScan->endMZ_ = filterLine.scanRangeMax_[ filterLine.scanRangeMax_.size() - 1 ];
	}
	else {
		// SRM: record range of q3 transition lists
		// start mz is average of first transition range
		// end mz is average of last transition range
		curScan->startMZ_ = (filterLine.transitionRangeMin_[0] + 
							 filterLine.transitionRangeMin_ [filterLine.transitionRangeMin_.size() - 1 ])
							 / 2;
		curScan->endMZ_ =  (filterLine.transitionRangeMax_[0] + 
							 filterLine.transitionRangeMax_ [filterLine.transitionRangeMax_.size() - 1 ])
							 / 2;
	}
	// get additional header information through Xcalibur:
	// retention time
	// min/max observed mz
	// total ion current
	// base peak mz and intensity
	// precursor (JMT: if avalible from interface, else use filter line info)

	long numDataPoints = -1; // points in both the m/z and intensity arrays
	double retentionTimeInMinutes = -1;
	long channel; // unused
	long uniformTime; // unused
	double frequency; // unused

	xrawfile2_->GetScanHeaderInfoForScanNum(
		scannumber, 
		&numDataPoints, 
		&retentionTimeInMinutes, 
		&(curScan->minObservedMZ_),
		&(curScan->maxObservedMZ_),
		&(curScan->totalIonCurrent_),
		&(curScan->basePeakMZ_),
		&(curScan->basePeakIntensity_),
		&channel, // unused
		&uniformTime, // unused
		&frequency // unused
		);

	// NOTE! the returned numDataPoints is invalid!
	// use the value from GetMassListFromScanNum below

	// record the retention time
	curScan->retentionTimeInSec_ = retentionTimeInMinutes * 60.0;

	// if ms level 2 or above, get precursor info
	if (curScan->msLevel_ > 1)  {
		getPrecursorInfo(*curScan, scannumber, filterLine);
	}


	//
	// get the m/z intensity pairs list for the current scan
	// !and correct min/max observed m/z info here!
	//

	
	curScan->minObservedMZ_ = 0;
	curScan->maxObservedMZ_ = 0;

	if (numDataPoints != 0) {


		// cout << "reading data points for scan " << curScanNum_ << endl;
		VARIANT varMassList;
		// initiallize variant to VT_EMPTY
		VariantInit(&varMassList);

		VARIANT varPeakFlags; // unused
		// initiallize variant to VT_EMPTY
		VariantInit(&varPeakFlags);

		// set up the parameters to read the scan
		// TODO make centroid parameter user customizable
		long dataPoints = 0;
		long scanNum = scannumber;
		LPCTSTR szFilter = NULL;		// No filter
		long intensityCutoffType = 0;		// No cutoff
		long intensityCutoffValue = 0;	// No cutoff
		long maxNumberOfPeaks = 0;		// 0 : return all data peaks
		double centroidPeakWidth = 0;		// No centroiding



		// record centroiding info
		//
		// scan may have been centroided at accquision time,
		// rather than conversion time (now)
		// even if user didn't request it.
		if ( (doCentroid_ /* && curScan->msLevel_ > 1 */)
			||
			filterLine.scanData_ == CENTROID) {
				curScan->isCentroided_ = true;	
		}

		// Note: special case for FT centroiding, contributed from Matt Chambers
		if (doCentroid_ /* && (curScan->msLevel_ > 1) */ && (curScan->analyzer_ == FTMS)) {
			// use GetLabelData to workaround bug in Thermo centroiding of FT profile data

			VARIANT varLabels;
			_variant_t vSpecData;
			xrawfile2_->GetLabelData(&varLabels, NULL, &scanNum);
			vSpecData.Attach(varLabels);
			dataPoints = vSpecData.parray->rgsabound[0].cElements;

			// record the number of data point (allocates memory for arrays)
			curScan->setNumDataPoints(dataPoints);

			double* pdval = (double*) vSpecData.parray->pvData;
			for(long i=0; i < dataPoints; ++i) {
				curScan->mzArray_[i] = (double) pdval[(i*6)+0];
				curScan->intensityArray_[i] = (double) pdval[(i*6)+1];
				//added by zhangjy 
				curScan->baseline_[i]=(double) pdval[(i*6)+3];
				curScan->noise_[i]=(double) pdval[(i*6)+4];
				//on 2011.1.5
			}
		} else {
			// get peaks with GetMassListFromScanNum as usual

			// do centroid should have no effect on already centroided data
			// but make sure to turn off centroiding for ms1 scans!
			bool centroidThisScan = doCentroid_;
			/*
			if (curScan->msLevel_ < 2) {
			centroidThisScan = false;
			}
			*/

			xrawfile2_->GetMassListFromScanNum(
				&scanNum,
				szFilter,			 // filter
				intensityCutoffType, // intensityCutoffType
				intensityCutoffValue, // intensityCutoffValue
				maxNumberOfPeaks,	 // maxNumberOfPeaks
				centroidThisScan,		// centroid result?
				&centroidPeakWidth,	// centroidingPeakWidth
				&varMassList,		// massList
				&varPeakFlags,		// peakFlags
				&dataPoints);		// array size

			// record the number of data point (allocates memory for arrays)
			curScan->setNumDataPoints(dataPoints);
			/*
			if (dataPoints != curScan->getNumDataPoints()) {
			cerr << "exiting with error at " 
			<< __FILE__ << ", line "
			<< __LINE__ << endl;
			exit(-1);
			}
			*/
			// Get a pointer to the SafeArray
			SAFEARRAY FAR* psa = varMassList.parray;
			DataPeak* pDataPeaks = NULL;
			SafeArrayAccessData(psa, (void**)(&pDataPeaks));

			// record mass list information in scan object
			if( doCentroid_ )
			{ // If we centroided the data we need to correct the basePeak m/z and intensity 
			  //(since GetScanHeaderInfoForScanNum returns information relevant to the profile data)!!
				curScan->basePeakIntensity_ = 0;
				for (long j=0; j<dataPoints; j++) {
					curScan->mzArray_[j] = pDataPeaks[j].dMass;
					curScan->intensityArray_[j] = pDataPeaks[j].dIntensity;
					if( pDataPeaks[j].dIntensity > curScan->basePeakIntensity_ )
					{
						curScan->basePeakMZ_ = pDataPeaks[j].dMass;
						curScan->basePeakIntensity_ = pDataPeaks[j].dIntensity;
					}
				}
			}
			else
			{
				for (long j=0; j<dataPoints; j++) {
					curScan->mzArray_[j] = pDataPeaks[j].dMass;
					curScan->intensityArray_[j] = pDataPeaks[j].dIntensity;
				}
			}

			// cleanup
			SafeArrayUnaccessData(psa); // Release the data handle
			VariantClear(&varMassList); // Delete all memory associated with the variant
			VariantClear(&varPeakFlags); // and reinitialize to VT_EMPTY

			if( varMassList.vt != VT_EMPTY ) {
				SAFEARRAY FAR* psa = varMassList.parray;
				varMassList.parray = NULL;
				SafeArrayDestroy( psa ); // Delete the SafeArray
			}

			if(varPeakFlags.vt != VT_EMPTY ) {
				SAFEARRAY FAR* psa = varPeakFlags.parray;
				varPeakFlags.parray = NULL;
				SafeArrayDestroy( psa ); // Delete the SafeArray
			}
		}

		// !!
		// Fix to overcome bug in ThermoFinnigan library GetScanHeaderInfoForScanNum() function
		// !!
		if (dataPoints > 0) {
			// don't do this on an empty scan!
			curScan->minObservedMZ_ = curScan->mzArray_[0];
			curScan->maxObservedMZ_ = curScan->mzArray_[dataPoints-1];
		}
	} // end 'not empty scan'
	else {
		// if empty scan:
		if (verbose_) 
		{
			//cout << "Note: empty scan detected (scan # " << scannumber << ")" << endl;
			printf("Note: empty scan detected (scan # %d\n", scannumber);
			//cout.flush();
			//flush();
		}
	}

	return curScan;
}
int ThermoInterface::ReselectPIons(DataPeak *PKL,int dpnum,double curpmz,int ch)
{
	if(ch<=0) return -1;	
	double dm=PRE_DM_10*curpmz/1e6;
	int idx=LocatePmz(PKL,dpnum,curpmz-dm,curpmz+dm,curpmz);
	if(idx==-1) return -1;
	int last1=-1;
	double maxInt=0;
	vector<int> isoIDX;	
	size_t i;
	int iso=1;
	double step=1.0033548/ch;
	double Lim=4*step+0.05;
	double cen=curpmz-step;
	double cL=cen-dm;
	double cH=cen+dm;
	for(i=idx;i>0;i--)
	{
		if(PKL[i].dMass<curpmz-Lim) break;
		if(PKL[i].dMass>cH) continue;
		else if(PKL[i].dMass<cL)
		{
			if(last1==-1) break;
			isoIDX.push_back(last1);
			cen-=step;
			cL=cen-dm;
			cH=cen+dm;
			last1=-1;
			maxInt=0;
			iso++;
		}
		else if(maxInt<PKL[i].dIntensity)
		{
			last1=i;
			maxInt=PKL[i].dIntensity;
		}
	}

	size_t sg=isoIDX.size();
	if(sg<=0) return idx;
	double ISOT[MAX_ISO];
	GetIsoDis(curpmz*ch,ISOT);
	for(i=0;i<sg;i++)
	{
		double ratioT=ISOT[0]/ISOT[i+1];
		double ratioE=PKL[isoIDX[i]].dIntensity/PKL[idx].dIntensity;
		if(ratioE>1.5*ratioT||ratioE<0.5*ratioT) break;
	}

	if(i==0) return idx;
	return isoIDX[i-1];
}

int ThermoInterface::TFindL(DataPeak *PKL,int b,int e,double fmz)
{
	if(e<b) return -1;
	double dm=PRE_DM_10*fmz/1e6;
	if(PKL[e].dMass<fmz-dm) return -1;
	if(PKL[b].dMass>fmz+dm) return -1;
	if(e-b<=1)
	{
		//if(PKL[b].dMass>fmz+dm||PKL[b].dMass<fmz-dm) return -1;
		return b;	
	}
	int idx=(b+e)/2;	
	if(PKL[idx].dMass>fmz)
	{
		return TFindL(PKL,b,idx,fmz);
	}
	else return TFindL(PKL,idx,e,fmz);
}

int ThermoInterface::TFindH(DataPeak *PKL,int b,int e,double fmz)
{
	if(e<b) return -1;
	double dm=PRE_DM_10*fmz/1e6;
	if(PKL[e].dMass<fmz-dm) return -1;
	if(PKL[b].dMass>fmz+dm) return -1;
	if(e-b<=1)
	{
		//if(PKL[e].dMass>fmz+dm||PKL[e].dMass<fmz-dm) return -1;
		return e;
	}
	int idx=(b+e)/2;
	if(PKL[idx].dMass>fmz)
	{
		return TFindH(PKL,b,idx,fmz);
	}
	else return TFindH(PKL,idx,e,fmz);
}

//Get the most high in the limit range
//lasted modifide on 2012.2.7,retrun to the original
int ThermoInterface::LocatePmz(DataPeak *PKL,int n,double begin,double end,double pmz)
{
	int i,middle=-1;	
	if(n<=0) return -1;
	double MaxInt=0;

	int iB=TFindL(PKL,0,n-1,begin);
	if(iB==-1) return -1;
	int iE=TFindH(PKL,0,n-1,end);
	if(iE==-1) return -1;
	for(i=iB;i<=iE;i++)
	{		
		if(PKL[i].dIntensity<=1e-6) continue;		
		if(PKL[i].dMass<begin) continue;
		if(PKL[i].dMass>end) break;
		if(MaxInt<PKL[i].dIntensity)
		{
			MaxInt=PKL[i].dIntensity;
			middle=i;
		}
	}
	return middle;
}

Scan* ThermoInterface::getScan(void) 
{
	if (!firstTime_) 
	{
		++curScanNum_;
		if (curScanNum_ > lastScanNumber_) 
		{
			// we're done
			return NULL;
		}
	} 
	else 
	{
		firstTime_ = false;
	}

	Scan* curScan = new Scan();
	curScan->isThermo_ = true;
	curScan->curScanNum=curScanNum_;

	//// test the "scan event" call
	//// gives a more through scan filter line
	//BSTR bstrScanEvent = NULL;
	//xrawfile2_->GetScanEventForScanNum(curScanNum_, &bstrScanEvent);
	//SysFreeString(bstrScanEvent);

	// Get the scan filter
	// (ex: "ITMS + c NSI Full ms [ 300.00-2000.00]")
	BSTR bstrFilter = NULL;
	xrawfile2_->GetFilterForScanNum(curScanNum_, &bstrFilter);
	curScan->thermoFilterLine_ = convertBstrToString(bstrFilter);
	SysFreeString(bstrFilter);
	FilterLine filterLine;
	if (!filterLine.parse(curScan->thermoFilterLine_)) 
	{
		//cerr << "error parsing filter line. exiting." << endl;
		//cerr << "line: " << curScan->thermoFilterLine_ << endl;
		printf("error parsing filter line. exiting.\n");
		printf("line: %d\n" ,curScan->thermoFilterLine_ );
		exit(-1);
	}

	// we should now have:
	// msLevel
	// polarity
	// scanType (full, zoom, etc)
	// ionization type
	// analyzer
	// scan range
	// parent mass and CID energy, if MS >=2 

	// record msLevel from filter line
	curScan->msLevel_ = filterLine.msLevel_;

	curScan->curScanNum=curScanNum_;//added by zhangjy on 2011.1.4

	// record polarity from filter line
	curScan->polarity_ = filterLine.polarity_;

	// record analyzer from filter line
	curScan->analyzer_ = filterLine.analyzer_;

	// record ionization from filter line
	curScan->ionization_ = filterLine.ionizationMode_;

	// record scan type from filter line
	// (zoom, full, srm, etc)
	curScan->scanType_ = filterLine.scanType_;

	// record activation (CID, etc)
	// check FilterLine: this may be default set to CID now
	if (filterLine.activationMethod_ != ACTIVATION_UNDEF) 
	{
		curScan->activation_ = filterLine.activationMethod_;
	}
	else 
	{
		curScan->activation_ = CID;
	}

	// record scan ranges from filter line
	// Note: SRM fills this with range of q3 transtion lists
	if (curScan->scanType_ != SRM) 
	{
		curScan->startMZ_ = filterLine.scanRangeMin_[ filterLine.scanRangeMin_.size() - 1 ];
		curScan->endMZ_ = filterLine.scanRangeMax_[ filterLine.scanRangeMax_.size() - 1 ];
	}
	else 
	{
		// SRM: record range of q3 transition lists
		// start mz is average of first transition range
		// end mz is average of last transition range
		curScan->startMZ_ = (filterLine.transitionRangeMin_[0] + 
							 filterLine.transitionRangeMin_ [filterLine.transitionRangeMin_.size() - 1 ])
							 / 2;
		curScan->endMZ_ =  (filterLine.transitionRangeMax_[0] + 
							 filterLine.transitionRangeMax_ [filterLine.transitionRangeMax_.size() - 1 ])
							 / 2;
	}
	// get additional header information through Xcalibur:
	// retention time
	// min/max observed mz
	// total ion current
	// base peak mz and intensity
	// precursor (JMT: if avalible from interface, else use filter line info)

	long numDataPoints = -1; // points in both the m/z and intensity arrays
	double retentionTimeInMinutes = -1;
	long channel; // unused
	long uniformTime; // unused
	double frequency; // unused

	xrawfile2_->GetScanHeaderInfoForScanNum(
		curScanNum_, 
		&numDataPoints, 
		&retentionTimeInMinutes, 
		&(curScan->minObservedMZ_),
		&(curScan->maxObservedMZ_),
		&(curScan->totalIonCurrent_),
		&(curScan->basePeakMZ_),
		&(curScan->basePeakIntensity_),
		&channel, // unused
		&uniformTime, // unused
		&frequency // unused
		);

	// NOTE! the returned numDataPoints is invalid!
	// use the value from GetMassListFromScanNum below

	// record the retention time
	curScan->retentionTimeInSec_ = retentionTimeInMinutes * 60.0;

	// if ms level 2 or above, get precursor info
	if (curScan->msLevel_ > 1) 
	{
		getPrecursorInfo(*curScan, curScanNum_, filterLine);
	}
	else//added by zhangjy on 2011.6.10
	{
		getStatusData(curScan);
		getTailerData(curScan);
	}

	//
	// get the m/z intensity pairs list for the current scan
	// !and correct min/max observed m/z info here!
	//

	
	curScan->minObservedMZ_ = 0;
	curScan->maxObservedMZ_ = 0;

	if (numDataPoints != 0) 
	{
		// cout << "reading data points for scan " << curScanNum_ << endl;
		VARIANT varMassList;
		// initiallize variant to VT_EMPTY
		VariantInit(&varMassList);

		VARIANT varPeakFlags; // unused
		// initiallize variant to VT_EMPTY
		VariantInit(&varPeakFlags);

		// set up the parameters to read the scan
		// TODO make centroid parameter user customizable
		long dataPoints = 0;
		long scanNum = curScanNum_;
		LPCTSTR szFilter = NULL;		// No filter
		long intensityCutoffType = 0;		// No cutoff
		long intensityCutoffValue = 0;	// No cutoff
		long maxNumberOfPeaks = 0;		// 0 : return all data peaks
		double centroidPeakWidth = 0;		// No centroiding

		// record centroiding info
		//
		// scan may have been centroided at accquision time,
		// rather than conversion time (now)
		// even if user didn't request it.
		if ( (doCentroid_ /* && curScan->msLevel_ > 1 */)
			||
			filterLine.scanData_ == CENTROID) 
		{
				curScan->isCentroided_ = true;	
		}

		// Note: special case for FT centroiding, contributed from Matt Chambers
		if (doCentroid_ /* && (curScan->msLevel_ > 1) */ && (curScan->analyzer_ == FTMS)) 
		{
			// use GetLabelData to workaround bug in Thermo centroiding of FT profile data

			VARIANT varLabels;
			_variant_t vSpecData;
			xrawfile2_->GetLabelData(&varLabels, NULL, &scanNum);
			vSpecData.Attach(varLabels);
			dataPoints = vSpecData.parray->rgsabound[0].cElements;

			// record the number of data point (allocates memory for arrays)
			curScan->setNumDataPoints(dataPoints);

			double* pdval = (double*) vSpecData.parray->pvData;
			for(long i=0; i < dataPoints; ++i)
			{
				curScan->mzArray_[i] = (double) pdval[(i*6)+0];
				curScan->intensityArray_[i] = (double) pdval[(i*6)+1];
				//added by zhangjy 
				curScan->baseline_[i]=(double) pdval[(i*6)+3];
				curScan->noise_[i]=(double) pdval[(i*6)+4];
				//on 2011.1.5
			}
		} 
		else 
		{
			// get peaks with GetMassListFromScanNum as usual

			// do centroid should have no effect on already centroided data
			// but make sure to turn off centroiding for ms1 scans!
			bool centroidThisScan = doCentroid_;
			/*
			if (curScan->msLevel_ < 2) {
			centroidThisScan = false;
			}
			*/

			xrawfile2_->GetMassListFromScanNum(
				&scanNum,
				szFilter,			 // filter
				intensityCutoffType, // intensityCutoffType
				intensityCutoffValue, // intensityCutoffValue
				maxNumberOfPeaks,	 // maxNumberOfPeaks
				centroidThisScan,		// centroid result?
				&centroidPeakWidth,	// centroidingPeakWidth
				&varMassList,		// massList
				&varPeakFlags,		// peakFlags
				&dataPoints);		// array size

			// record the number of data point (allocates memory for arrays)
			curScan->setNumDataPoints(dataPoints);
			/*
			if (dataPoints != curScan->getNumDataPoints()) {
			cerr << "exiting with error at " 
			<< __FILE__ << ", line "
			<< __LINE__ << endl;
			exit(-1);
			}
			*/
			// Get a pointer to the SafeArray
			SAFEARRAY FAR* psa = varMassList.parray;
			DataPeak* pDataPeaks = NULL;
			SafeArrayAccessData(psa, (void**)(&pDataPeaks));

			// record mass list information in scan object
			if( doCentroid_ )
			{ // If we centroided the data we need to correct the basePeak m/z and intensity 
			  //(since GetScanHeaderInfoForScanNum returns information relevant to the profile data)!!
				curScan->basePeakIntensity_ = 0;
				for (long j=0; j<dataPoints; j++) 
				{
					curScan->mzArray_[j] = pDataPeaks[j].dMass;
					curScan->intensityArray_[j] = pDataPeaks[j].dIntensity;
					if( pDataPeaks[j].dIntensity > curScan->basePeakIntensity_ )
					{
						curScan->basePeakMZ_ = pDataPeaks[j].dMass;
						curScan->basePeakIntensity_ = pDataPeaks[j].dIntensity;
					}
				}
			}
			else
			{
				for (long j=0; j<dataPoints; j++) 
				{
					curScan->mzArray_[j] = pDataPeaks[j].dMass;
					curScan->intensityArray_[j] = pDataPeaks[j].dIntensity;
				}
			}

				// cleanup
			SafeArrayUnaccessData(psa); // Release the data handle
			VariantClear(&varMassList); // Delete all memory associated with the variant
			VariantClear(&varPeakFlags); // and reinitialize to VT_EMPTY

			if( varMassList.vt != VT_EMPTY ) 
			{
					SAFEARRAY FAR* psa = varMassList.parray;
					varMassList.parray = NULL;
					SafeArrayDestroy( psa ); // Delete the SafeArray
			}

			if(varPeakFlags.vt != VT_EMPTY )
			{
					SAFEARRAY FAR* psa = varPeakFlags.parray;
					varPeakFlags.parray = NULL;
					SafeArrayDestroy( psa ); // Delete the SafeArray
			}
	    }

		// !!
		// Fix to overcome bug in ThermoFinnigan library GetScanHeaderInfoForScanNum() function
		// !!
		if (dataPoints > 0) 
		{
			// don't do this on an empty scan!
			curScan->minObservedMZ_ = curScan->mzArray_[0];
			curScan->maxObservedMZ_ = curScan->mzArray_[dataPoints-1];
		}
	} // end 'not empty scan'
	else 
	{
		// if empty scan:
		if (verbose_) 
		{
			//cout << "Note: empty scan detected (scan # " << curScanNum_ << ")" << endl;
			printf("Note: empty scan detected (scan # %d\n",curScanNum_);
			//cout.flush();
			//flush();
		}
	}
	//在返回之前重新进行处理
	curScan->IsCalibrate=false;	
	if(IsCalibrate)
	{	
		if(curScan->msLevel_==1)
		{
			if(MType==MODEL_SVM) CalibeateMS_SVM(curScan);
			else if(MType==MODEL_LIN) CalibeateMS_LIN(curScan);
		}
		else if(curScan->msLevel_==2)
		{
			long idx=seekScan(curScan->curScanNum,0,totalpmzscan);
			if(idx!=-1)
			{
				curScan->precursorMZ_=CalPmz[idx];
			}
		}
	}
	return curScan;
}

//fea[0].value=Isomz[0];
//	fea[1].value=RT;
//	fea[2].value=tic;
//	fea[3].value=log10(IsotopicE[0]+1e-4);
//	fea[4].value=sqrt(S_N[0]);	
//	fea[5].value=0;
//	fea[6].value=ch;
//	fea[7].value=IsoNum;
//	size_t sg=ParSel.size();
//	for(size_t i=0;i<sg;i++)
//	{
//		fea[8+i].value=FTStatus[ParSel[i]];
//	}	

void ThermoInterface::CalibeateMS_SVM(Scan *curScan)
{
	size_t i,VALFea=statusIDX.size();		
	int FEA_NUM=FIXED_PAR_NUM_SVM+VALFea;
	if(trailerIDX[0]>0) FEA_NUM++;
	if(trailerIDX[1]>0) FEA_NUM++;
	double *fea;
	fea=new double[FEA_NUM];	
	fea[1]=curScan->retentionTimeInSec_/60;	
	size_t pnum=curScan->getNumDataPoints();
	fea[2]=0;
	for(i=0;i<pnum;i++) fea[2]+=curScan->intensityArray_[i]/TIC_SCALE;		
	size_t j=FIXED_PAR_NUM_SVM;
	if(trailerIDX[0]>0)
	{
		fea[j]=atof(curScan->tailer_par_value[trailerIDX[0]].c_str());
		j++;
	}
	if(trailerIDX[1]>0)
	{
		fea[j]=atof(curScan->tailer_par_value[trailerIDX[1]].c_str());
		j++;
	}

	for(i=0;i<VALFea;i++)
	{			
		fea[j]=atof(curScan->status_par_value[statusIDX[i]].c_str());
		j++;
	}		
	for(i=0;i<pnum;i++)
	{
		fea[0]=curScan->mzArray_[i];
		fea[3]=log10(curScan->intensityArray_[i]+1e-4);
		fea[5]=curScan->intensityArray_[i]/curScan->basePeakIntensity_;		
		fea[4]=curScan->noise_[i]+curScan->baseline_[i];
		if(fea[4]>1e-2)fea[4]=sqrt(curScan->intensityArray_[i]/fea[4]);	
		curScan->mzArray_[i]=CalModel->Predict(fea,FEA_NUM);		
	}
	curScan->IsCalibrate=true;	
	delete []fea;
}


void ThermoInterface::CalibeateMS_LIN(Scan *curScan)
{
	size_t i,VALFea=statusIDX.size();		
	int FEA_NUM=FIXED_PAR_NUM_LIN+VALFea;
	if(trailerIDX[0]>0) FEA_NUM++;
	if(trailerIDX[1]>0) FEA_NUM++;
	double *fea;
	fea=new double[FEA_NUM];	
	fea[1]=curScan->retentionTimeInSec_/60;	
	size_t pnum=curScan->getNumDataPoints();
	fea[2]=0;
	for(i=0;i<pnum;i++) fea[2]+=curScan->intensityArray_[i]/TIC_SCALE;			
	size_t j=FIXED_PAR_NUM_LIN;
	if(trailerIDX[0]>0)
	{
		fea[j]=atof(curScan->tailer_par_value[trailerIDX[0]].c_str());
		j++;
	}
	if(trailerIDX[1]>0)
	{
		fea[j]=atof(curScan->tailer_par_value[trailerIDX[1]].c_str());
		j++;
	}

	for(i=0;i<VALFea;i++)
	{			
		fea[j]=atof(curScan->status_par_value[statusIDX[i]].c_str());
		j++;
	}		
	for(i=0;i<pnum;i++)
	{
		fea[0]=curScan->mzArray_[i];
		fea[3]=log10(curScan->intensityArray_[i]+1e-4);
		fea[4]=curScan->noise_[i]+curScan->baseline_[i];
		fea[5]=curScan->intensityArray_[i]/curScan->basePeakIntensity_;//iso	
		if(fea[4]>1e-2)	fea[4]=sqrt(curScan->intensityArray_[i]/fea[4]);	
		curScan->mzArray_[i]=LM.Predict(fea,FEA_NUM);		
	}
	curScan->IsCalibrate=true;	
	delete []fea;
}

double ThermoInterface::GetSampleDm(double mass)
{

	if(instrumentInfo_.instrumentModel_==LTQ_FT||instrumentInfo_.instrumentModel_==LTQ_FT_ULTRA) return 6.06e-9*mass*mass-7.063e-11*mass+4.085e-8;//for FT
	else if(instrumentInfo_.instrumentModel_==LTQ_ORBITRAP||instrumentInfo_.instrumentModel_==LTQ_ORBITRAP_DISCOVERY||
		instrumentInfo_.instrumentModel_==LTQ_ORBITRAP_XL||instrumentInfo_.instrumentModel_==LTQ_ORBITRAP_VELOS)
		return 1.12e-009*mass*mass+2.204e-006*mass-0.0003366;//for Orbitrap
	else return 0;
}

void ThermoInterface::SortMerge(vector<double> &tmz,vector<double> &tint)
{
	//sort at first
	size_t i,k,sg=tmz.size();
	for(i=0;i<sg;i++)
	{
		for(k=i+1;k<sg;k++)
		{
			if(tmz[i]>tmz[k]) 
			{
				double mt=tmz[i];
				tmz[i]=tmz[k];
				tmz[k]=mt;

				mt=tint[i];
				tint[i]=tint[k];
				tint[k]=mt;
			}
		}
	}
	//merge the neibour peaks
	k=0;
	vector<double> Merg_mz;
	vector<double> Merg_dint;
	while(k<sg)
	{	
		double dm_lim=GetSampleDm(tmz[k])*5;//half peaks width
		for(i=k+1;i<sg;i++)
		{
			double dm=tmz[i]-tmz[k];
			if(dm>dm_lim) break;
		}
		double dMass=0;
		double dIntensity=0;
		for(;k<i;k++)
		{
			dMass+=tmz[k]*tint[k];
			dIntensity+=tint[k];
		}
		if(dIntensity>1e-2)dMass/=dIntensity;
		Merg_mz.push_back(dMass);	
		Merg_dint.push_back(dIntensity);
	}
	tmz.clear();
	tint.clear();
	sg=Merg_mz.size();
	for(i=0;i<sg;i++)
	{
		tmz.push_back(Merg_mz[i]);
		tint.push_back(Merg_dint[i]);
	}
}

// get precursor m/z, collision energy, precursor charge, and precursor intensity
bool ThermoInterface::getPrecursorInfo(Scan& scan, long scanNumber, FilterLine& filterLine) 
{
	// TODO: assert scan->msLevel_ > 1

	if (scanNumber == 1) 
	{
		// if this is the first scan, only use the info from the filter lin
		// (API calls are no use, as there's no precursor scan)
		// An example of a first scan WITH filterline precursor info:
		// "+ c NSI sid=5.00  SRM ms2 580.310@cid25.00 [523.955-523.965, 674.405-674.415, 773.295-773.305]"

		scan.accuratePrecursorMZ_ = false;
		// use the low-precision parent mass in the filter line
		inaccurateMasses_++;
		scan.precursorMZ_ = filterLine.cidParentMass_[filterLine.cidParentMass_.size() - 1];
		scan.accuratePrecursorMZ_ = false;		

		// use the low-precision collision energy recorded in the filter line
		scan.collisionEnergy_ = filterLine.cidEnergy_[filterLine.cidEnergy_.size() - 1];

		scan.precursorCharge_ = -1; // undetermined
		chargeCounts_[0]++; // with "charge 0" being undetermined

		scan.precursorIntensity_ = 0; // undetermined

		return true;
	}



	// first try with the more direct GetPrecursorInfoFromScanNum call;
	// if this fails, continue with older methods
	if (IXRawfileVersion_ > 2)
	{
		try
		{
			VARIANT vPrecursorInfos;
			VariantInit(&vPrecursorInfos);
			long nPrecursorInfos = 0;


			// Get the precursor scan information
			XRAWFILE2Lib::IXRawfile3 * ix3ptr = (XRAWFILE2Lib::IXRawfile3 *)xrawfile2_.GetInterfacePtr();

			ix3ptr->GetPrecursorInfoFromScanNum(scanNumber,
				&vPrecursorInfos, 
				&nPrecursorInfos);


			// for now, exit if anything more than 1 precursor found
			if (nPrecursorInfos > 1) 
			{
				//cerr << "halting; multiple precursors detected.  please contact developer." << endl;
				printf("halting; multiple precursors detected.  please contact developer.\n");
				exit(1);
			}

			// Access the safearray buffer
			BYTE* pData;
			SafeArrayAccessData(vPrecursorInfos.parray, (void**)&pData);

			for (int i=0; i < nPrecursorInfos; ++i) {
				// Copy the scan information from the safearray buffer
				//PrecursorInfo info;
				XRAWFILE2Lib::MS_PrecursorInfo precursorInfo;
				memcpy(&precursorInfo, 
					pData + i * sizeof(XRAWFILE2Lib::MS_PrecursorInfo), 
					sizeof(XRAWFILE2Lib::MS_PrecursorInfo));

				// set the scan info
				scan.precursorCharge_ = precursorInfo.nChargeState;
				scan.precursorScanNumber_ = precursorInfo.nScanNumber;
				scan.precursorMZ_ = precursorInfo.dMonoIsoMZ;
				if (precursorInfo.dMonoIsoMZ != 0) {
					scan.accuratePrecursorMZ_ = true;
					accurateMasses_++;
					getPreInfoCount_++;
				}
			}

			SafeArrayUnaccessData(vPrecursorInfos.parray);
		}
		catch (...)
		{
			//cerr << "There was a problem while getting the precursor scan information (scan " <<  
			//	scanNumber <<"); exiting" << endl;
			//cerr << "please contact developer with this error message." << endl;
			printf("There was a problem while getting the precursor scan information (scan ");
			printf("%d ); exiting\n",scanNumber);
			printf("please contact developer with this error message.\n");
		}
	}

	// hopefully we have everything except precursor intensity and collision energy here.
	// if not, fall back on the older versions.


	VARIANT varValue;
	VariantInit(&varValue);

	//added by zhangjy 
	xrawfile2_->GetTrailerExtraValueForScanNum(curScanNum_, "MS2 Isolation Width:" , &varValue);
	if( varValue.vt == VT_R4 ) 
	{
		scan.IsolationWidth_= varValue.fltVal;
	}
	else if( varValue.vt == VT_R8 ) 
	{
		scan.IsolationWidth_ = varValue.dblVal;
	}
	else if ( varValue.vt != VT_ERROR ) 
	{
		//cerr << "Scan: " << curScanNum_ << " MS level: " << scan.msLevel_ 
		//	<< " unexpected type when looking for precursorMz\n";
		//exit(-1);

		return false;
	}

	/*
	// compare new vs old API
	double oldMZ = 0;
	xrawfile2_->GetTrailerExtraValueForScanNum(curScanNum_, "Monoisotopic M/Z:" , &varValue);

	if( varValue.vt == VT_R4 ) {
		oldMZ = varValue.fltVal;
	}
	else if( varValue.vt == VT_R8 ) {
		oldMZ = varValue.dblVal;
	}
	else if ( varValue.vt != VT_ERROR ) {
		cerr << "Scan: " << curScanNum_ << " MS level: " << scan.msLevel_ 
			<< " unexpected type when looking for precursorMz\n";
		exit(-1);
	}
	if (oldMZ != scan.precursorMZ_) {
		cerr << "scan: " << setprecision(11) << curScanNum_ << " old mz=" << oldMZ << " new mz=" << scan.precursorMZ_ << endl;
	}
	*/

	//
	// precursor m/z fallbacks
	// try the Thermo API first, then falling back on the value extracted from the filter line
	//

	if (scan.precursorMZ_ <= 0) { // restore if getPrec() call works in the future
//	if (1) {
		// we'll try to get the "accurate mass" recorded by the machine, and fall back to the filter line value otherwise.
		scan.accuratePrecursorMZ_ = false;

		// see if we can get an accurate value from the machine

		double precursorMZ = 0;


		// don't try to get Monoisotopic M/Z if the user has selected "force precursor from filter line"
		if (!forcePrecursorFromFilter_) {
			// ignore return value from this call
			xrawfile2_->GetTrailerExtraValueForScanNum(curScanNum_, "Monoisotopic M/Z:" , &varValue);

			if( varValue.vt == VT_R4 ) {
				precursorMZ = varValue.fltVal;
			}
			else if( varValue.vt == VT_R8 ) {
				precursorMZ = varValue.dblVal;
			}
			else if ( varValue.vt != VT_ERROR )
			{
			/*	cerr << "Scan: " << curScanNum_ << " MS level: " << scan.msLevel_ 
					<< " unexpected type when looking for precursorMz\n";*/
				//exit(-1);
				return false;
			}
		}

		scan.precursorMZ_ = filterLine.cidParentMass_[filterLine.cidParentMass_.size() - 1];
		if (precursorMZ > 0 && fabs(precursorMZ-scan.precursorMZ_)<=10.0) { // (note: this could only true if we tried and had sucess
			// with the monoisotopic call above.)
			// Sanity check to make sure mono mass is in ballpark of filter line mass.

			// sucess-- higher accuracy m/z recorded through API
			accurateMasses_++;
			scan.precursorMZ_ = precursorMZ;
			scan.accuratePrecursorMZ_ = true;
			oldAPICount_++;
		}
		else { 
			// use the low-precision parent mass in the filter line
			inaccurateMasses_++;
			//scan.precursorMZ_ = filterLine.cidParentMass_[filterLine.cidParentMass_.size() - 1];
			scan.accuratePrecursorMZ_ = false;
			filterLineCount_++;
		}
	}


	//
	// collision energy, trying the Thermo API first, then falling back on the value extracted from the filter line
	//

	double collisionEnergy = 0;
	VariantClear(&varValue);
	xrawfile2_->GetTrailerExtraValueForScanNum(curScanNum_, "API Source CID Energy:" , &varValue);
	if( varValue.vt == VT_R4 ) {
		// VT_R4: OLE float type
		collisionEnergy = varValue.fltVal;
	}		
	else if ( varValue.vt != VT_ERROR ) 
	{
		/*cerr << "Unexpected type when looking for CE\n";*/
		//exit(-1);
		return false;
	}

	if (collisionEnergy != 0) {
		// sucess-- collision energy recorded through API
		scan.collisionEnergy_ = collisionEnergy;
	}
	else {
		// use the low-precision collision energy recorded in the filter line
		scan.collisionEnergy_ = filterLine.cidEnergy_[filterLine.cidEnergy_.size() - 1];
	}

/*
	// test charge anyhow

	{
		VariantClear(&varValue);
		int precursorCharge = 0;
		xrawfile2_->GetTrailerExtraValueForScanNum(curScanNum_, "Charge State:" , &varValue);
		// First try to use the OCX
		if( varValue.vt == VT_I2 ) {
			precursorCharge = varValue.iVal;
		}

		if (precursorCharge != 0) {
			// sucess-- precursor charge recorded through API
			scan.precursorCharge_ = precursorCharge;
		}
		else {
			// no luck
			scan.precursorCharge_ = -1; // undetermined
		}
		if (precursorCharge == 0) {
			if (scan.precursorCharge_ != 0) {
				cerr << "old API was zero, but new call was " << scan.precursorCharge_ << endl;
			}
		}
		else {
			if (precursorCharge != scan.precursorCharge_) {
				cerr << "charge diff: old: " << precursorCharge << " new : " << scan.precursorCharge_ << endl;
			}
			else {
				cerr << "old/new agreed" << endl;
			}
		}
	}
*/

	int trailerPrecursorCharge = 0;
	if (scan.precursorCharge_ < 1) {
		// precursor charge state fallbacks, again trying first from the API, then resorting to filter line value
		VariantClear(&varValue);
		
		xrawfile2_->GetTrailerExtraValueForScanNum(curScanNum_, "Charge State:" , &varValue);
		// First try to use the OCX
		if( varValue.vt == VT_I2 ) {
			trailerPrecursorCharge = varValue.iVal;
		}

		if (trailerPrecursorCharge != 0) {
			// sucess-- precursor charge recorded through API
			scan.precursorCharge_ = trailerPrecursorCharge;
		}
		else {
			// no luck
			scan.precursorCharge_ = -1; // undetermined
		}
	}

	// track the counts of precursor charges
	// with "0" being undetermined
	if (scan.precursorCharge_ < 1) {
		chargeCounts_[0]++;
	}
	else {
		if ( (scan.precursorCharge_ + 1) > int(chargeCounts_.size())) {
			chargeCounts_.resize(scan.precursorCharge_ + 1, 0);
			//if (verbose_) {
			//	cout << "new max observed charge: " << scan.precursorCharge_ << "(scan " << curScanNum_ << ")" << endl;
			//	cout.flush();
			//}
		}
		chargeCounts_[scan.precursorCharge_]++;
	}


	//
	// precursor intensity determiniation
	//

	// go to precursor scan and try to find highest intensity around precursor m/z.
	// (could this be improved to handle more complex acquisition methods?)

	//if( numDataPoints != 0 ) { // if this isn't an empty scan
	VARIANT varMassList;
	VariantInit(&varMassList);	// initiallize variant to VT_EMPTY
	VARIANT varPeakFlags; // unused
	VariantInit(&varPeakFlags);	// initiallize variant to VT_EMPTY



	// set up the parameters to read the precursor scan
	long numPrecursorDataPoints = 0;
	LPCTSTR szFilter = "!d";		// First previous not-dependent scan
	long intensityCutoffType = 0;	// No cutoff
	long intensityCutoffValue = 0;	// No cutoff
	long maxNumberOfPeaks = 0;		// Return all data peaks
	BOOL centroidResult = FALSE;	// No centroiding of the precursor
	double centroidPeakWidth = 0;	// (see above: no centroiding)

	// cout << "reading scan " << curScanNum_ << endl;

	long curScanNum = scanNumber;

	// the goal is to get the parent scan's info
	xrawfile2_->GetPrevMassListFromScanNum(
		&curScanNum,
		szFilter,				// filter
		intensityCutoffType,	// intensityCutoffType
		intensityCutoffValue,	// intensityCutoffValue
		maxNumberOfPeaks,		// maxNumberOfPeaks
		centroidResult,			// centroidResult
		&centroidPeakWidth,		// centroidingPeakWidth
		&varMassList,			// massList
		&varPeakFlags,			// peakFlags -- unused
		&numPrecursorDataPoints);	// arraySize

	// record the precursor scan number
	scan.precursorScanNumber_ = curScanNum; // set during last xrawfile2 call

	// numPrecursorDataPoints: number of mass/intensity pairs for precursor scan
	// if empty, no precursor scan
	if (numPrecursorDataPoints != 0) {
		// Get a pointer to the SafeArray
		SAFEARRAY FAR* psa = varMassList.parray;
		DataPeak* pDataPeaks = NULL;
		// use our pDataPeaks array to access the SafeArray data.
		SafeArrayAccessData(psa, (void**)(&pDataPeaks) );
		double precursorIntensity = 0;

		// Find most intense peak around precursorMz:		
		// search the precursor's scan data near precursorMZ
		// to get the intensity.
		//for (long j=0; j<numPrecursorDataPoints; j++) {
		//	double dMass = pDataPeaks[j].dMass;

		//	if (fabs(scan.precursorMZ_ - dMass) < 0.05 ) {
		//		if( pDataPeaks[j].dIntensity > precursorIntensity ) {
		//			precursorIntensity = pDataPeaks[j].dIntensity;
		//		}
		//	}
		//	// stop if we've gone too far
		//	if (dMass - scan.precursorMZ_ > 0.05) {
		//		break;
		//	}
		//}
		//// record the precursor intensity
		//scan.precursorIntensity_ = precursorIntensity;
		int idx=ReselectPIons(pDataPeaks,numPrecursorDataPoints,scan.precursorMZ_,scan.precursorCharge_);
		if(idx!=-1) 
		{
			scan.precursorMZ_=pDataPeaks[idx].dMass;
			scan.precursorIntensity_=pDataPeaks[idx].dIntensity;
		}	
		else scan.precursorIntensity_=0;

		// Release the data handle
		SafeArrayUnaccessData(psa);
		// Delete all memory associated with the variant and
		// reinitialize to VT_EMPTY
		VariantClear(&varMassList);	
		VariantClear(&varPeakFlags);
	}

	// double check on deallocating memory
	if (varMassList.vt != VT_EMPTY) {
		SAFEARRAY FAR* psa = varMassList.parray;
		varMassList.parray = NULL;
		// Delete the SafeArray
		SafeArrayDestroy( psa );
	}

	if (varPeakFlags.vt != VT_EMPTY) {
		SAFEARRAY FAR* psa = varPeakFlags.parray;
		varPeakFlags.parray = NULL;
		// Delete the SafeArray
		SafeArrayDestroy(psa);
	}
	return true;
}

bool ThermoInterface::getStatusData(Scan *cs)
{
	long ScanNum=cs->curScanNum;
	double dStatusLogRT = 0.0;
	VARIANT varLabels;
	VariantInit(&varLabels);
	VARIANT varValues;
	VariantInit(&varValues);
	long nArraySize = 0;

	xrawfile2_->GetStatusLogForScanNum(ScanNum, 
	&dStatusLogRT, 
	&varLabels, 
	&varValues, 
	&nArraySize);
	if(nArraySize<=0) return false;

	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psaLabels = varLabels.parray;
	varLabels.parray = NULL;

	SAFEARRAY FAR* psaValues = varValues.parray;
	varValues.parray = NULL;

	BSTR* pbstrLabels = NULL;
	BSTR* pbstrValues = NULL;

	if(FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels))))
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		//printf("Failed to access labels array");
		return false;
	}

	if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		SafeArrayUnaccessData(psaValues);
		SafeArrayDestroy(psaValues);
		//printf("Failed to access values array");
		return false;
	}

	string sLabel;
	string sData;
	long i;
	for( i=0;i<nArraySize;i++ )
	{
		sLabel = convertBstrToString(pbstrLabels[i]);
		sData = convertBstrToString(pbstrValues[i]);			
	
		//printf("Name=%s\n",(LPCTSTR)sLabel);
		//printf("Value=%s\n",(LPCTSTR)sData);
		cs->status_par_name.push_back(sLabel);
		cs->status_par_value.push_back(sData);
	}
	// Delete the SafeArray
	SafeArrayUnaccessData(psaLabels);
	SafeArrayDestroy(psaLabels);
	SafeArrayUnaccessData(psaValues);
	SafeArrayDestroy(psaValues);	
	return true;
}

bool ThermoInterface::getTailerData(Scan *cs)
{
	long ScanNum=cs->curScanNum;
	double dStatusLogRT = 0.0;
	VARIANT varLabels;
	VariantInit(&varLabels);
	VARIANT varValues;
	VariantInit(&varValues);
	long nArraySize = 0;

	xrawfile2_->GetTrailerExtraForScanNum(ScanNum, 
	&varLabels, 
	&varValues, 
	&nArraySize);
	if(nArraySize<=0) return false;

	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psaLabels = varLabels.parray;
	varLabels.parray = NULL;

	SAFEARRAY FAR* psaValues = varValues.parray;
	varValues.parray = NULL;

	BSTR* pbstrLabels = NULL;
	BSTR* pbstrValues = NULL;

	if(FAILED(SafeArrayAccessData( psaLabels, (void**)(&pbstrLabels))))
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		//printf("Failed to access labels array");
		return false;
	}

	if( FAILED(SafeArrayAccessData( psaValues, (void**)(&pbstrValues) ) ) )
	{
		SafeArrayUnaccessData(psaLabels);
		SafeArrayDestroy(psaLabels);
		SafeArrayUnaccessData(psaValues);
		SafeArrayDestroy(psaValues);
		//printf("Failed to access values array");
		return false;
	}

	string sLabel;
	string sData;
	long i;
	for( i=0;i<nArraySize;i++ )
	{
		sLabel = convertBstrToString(pbstrLabels[i]);
		sData = convertBstrToString(pbstrValues[i]);			
	
		//printf("Name=%s\n",(LPCTSTR)sLabel);
		//printf("Value=%s\n",(LPCTSTR)sData);
		cs->tailer_par_name.push_back(sLabel);
		cs->tailer_par_value.push_back(sData);
	}
	// Delete the SafeArray
	SafeArrayUnaccessData(psaLabels);
	SafeArrayDestroy(psaLabels);
	SafeArrayUnaccessData(psaValues);
	SafeArrayDestroy(psaValues);	
	return true;
}

long ThermoInterface::getFirstScanNumber()
{
	return firstScanNumber_;
}

long ThermoInterface::GetLastScanNumber()
{
	return lastScanNumber_;
}

int ThermoInterface::MSLevel(long ScanNumber)
{
	if(ScanNumber<firstScanNumber_||ScanNumber>lastScanNumber_) return -1;
	BSTR bstrFilter = NULL;
	xrawfile2_->GetFilterForScanNum(ScanNumber, &bstrFilter);
	string thermoFilterLine = convertBstrToString(bstrFilter);
	SysFreeString(bstrFilter);
	FilterLine filterLine;
	if (!filterLine.parse(thermoFilterLine))  
	{
		//cerr << "error parsing filter line. exiting." << endl;
		//cerr << "line: " << thermoFilterLine << endl;	
		return -1;	
	}
	return filterLine.msLevel_;
}

int ThermoInterface::MSLevel()
{
	long testScan=curScanNum_;
	if(testScan<firstScanNumber_||testScan>lastScanNumber_) return -1;
	BSTR bstrFilter = NULL;
	xrawfile2_->GetFilterForScanNum(testScan, &bstrFilter);
	string thermoFilterLine = convertBstrToString(bstrFilter);
	SysFreeString(bstrFilter);
	FilterLine filterLine;
	if (!filterLine.parse(thermoFilterLine)) 
	{
	/*	cerr << "error parsing filter line. exiting." << endl;
		cerr << "line: " << thermoFilterLine << endl;*/
		return -1;
	}
	return filterLine.msLevel_;
}

bool ThermoInterface::SkipScan()
{
	curScanNum_++;
	if(curScanNum_>=lastScanNumber_) return false;
	return true;
}

//for each time , the getScan function will add 1 to curScanNum_
void ThermoInterface::SeekToScan(long targetscan)
{
	if(targetscan>lastScanNumber_||targetscan<firstScanNumber_)  return;
	else curScanNum_=targetscan-1;
}

#define DPNUM 32000
bool ThermoInterface::PreSuryMZINT()
{
	double x[DPNUM];
	double y[DPNUM];
	long ActNUM=0;
	bool centriodS=doCentroid_;
	doCentroid_=false;
	Scan *st;
	bool EndADD=false;
	for(long tscan=firstScanNumber_;tscan<lastScanNumber_;tscan++) 
	{
		if(MSLevel(tscan)!=1) continue;
		st=getScan(tscan);
		size_t as=st->getNumDataPoints();
		for(size_t i=1;i<as;i++)
		{
			double dm=st->mzArray_[i]-st->mzArray_[i-1];
			if(dm<0.05)
			{
				x[ActNUM]=st->mzArray_[i];
				y[ActNUM]=dm;
				ActNUM++;
				if(ActNUM>=DPNUM) 
				{
					EndADD=true;
					break;
				}
			}
		}
		delete st;
		if(EndADD) break;
	}
	MzIntModel[3]=gsl_two_step_reg(x,y,3,ActNUM,MzIntModel);
	doCentroid_=centriodS;
	return true;
}

void ThermoInterface::Convert(vector<size_t> &ColParIdx,vector<size_t> &SelParIdx,size_t tIDX[2])
{
	trailerIDX[0]=0;
	trailerIDX[1]=0;
	statusIDX.clear();	
	size_t sg=SelParIdx.size();

	for(size_t i=0;i<sg;i++)
	{
		if(SelParIdx[i]>1)
		{
			statusIDX.push_back(ColParIdx[SelParIdx[i]-2]);
		}
		else if(SelParIdx[i]==0) trailerIDX[0]=tIDX[0];
		else if(SelParIdx[i]==1) trailerIDX[1]=tIDX[1];		
	}
}

long ThermoInterface::seekScan(long scannum,size_t a,size_t b)
{
	if(b<a) return -1;
	if(b-a<=1) 
	{
		if(scanmap.at(a)==scannum)	return a;	
		else if(scanmap.at(b)==scannum) return b;
		else return -1;
	}
	int idx=(b+a)/2;
	if(scanmap.at(idx)==scannum) return idx;
	else if(scanmap.at(idx)>scannum)
	{
		return seekScan(scannum,a,idx);
	}
	else return seekScan(scannum,idx,b);
}

int ThermoInterface::countMS2()
{
	BSTR bstrFilter = NULL;
	string thermofilterline;
	FilterLine filterLine;
	int count=0;
	for(long scan=firstScanNumber_;scan<=lastScanNumber_;scan++)
	{
		xrawfile2_->GetFilterForScanNum(scan, &bstrFilter);
		thermofilterline= convertBstrToString(bstrFilter);
		SysFreeString(bstrFilter);
		 bstrFilter = NULL;
		if(filterLine.parse(thermofilterline)) 
		{
			if(filterLine.msLevel_==2) count++;		
		}
	}
	return count;
}

int ThermoInterface::countMS1()
{
	BSTR bstrFilter = NULL;
	string thermofilterline;
	FilterLine filterLine;
	int count=0;
	for(long scan=firstScanNumber_;scan<=lastScanNumber_;scan++)
	{
		xrawfile2_->GetFilterForScanNum(scan, &bstrFilter);
		thermofilterline= convertBstrToString(bstrFilter);
		SysFreeString(bstrFilter);
		bstrFilter = NULL;
		if(filterLine.parse(thermofilterline)) 
		{
			if(filterLine.msLevel_==1) count++;		
		}
	}
	return count;
}

void ThermoInterface::countMS(int CT[2])
{
	CT[0]=0;
	CT[1]=0;
	BSTR bstrFilter = NULL;
	string thermofilterline;
	FilterLine filterLine;
	for(long scan=firstScanNumber_;scan<=lastScanNumber_;scan++)
	{
		xrawfile2_->GetFilterForScanNum(scan,&bstrFilter);
		thermofilterline= convertBstrToString(bstrFilter);
		SysFreeString(bstrFilter);
		bstrFilter = NULL;
		if(filterLine.parse(thermofilterline)) 
		{
			if(filterLine.msLevel_==1) CT[0]++;	
			else if(filterLine.msLevel_==2) CT[1]++;	
		}
	}

}