// -*- mode: c++ -*-


/*
    File: Scan.cpp
    Description: instrument-independent scan representation.
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


#include "Scan.h"
#include "math.h"
#include <string>
#include <iostream>
#include <vector>

#include "GaussFit.h"
#define SQRTPI 1.7724538509055160272981674833411

using namespace std;

// From Spectrast, for Henry's centroiding:
// a peak - note that we use float's for intensities to save memory
typedef struct _peak {
	double mz;
	double intensity;
	string annotation;	
	string info;
} Peak;


void Scan::setNumDataPoints(int numDataPoints) {
	if (numDataPoints == 0) {
		numDataPoints_ = numDataPoints;
		return;
	}

	delete [] mzArray_;
	delete [] intensityArray_;

	delete []baseline_;
	delete []noise_;

	numDataPoints_ = numDataPoints;
	mzArray_ = new double[numDataPoints_];
	intensityArray_ = new double[numDataPoints_];
	baseline_=new double[numDataPoints_];
	noise_=new double[numDataPoints_];
}


void Scan::resetNumDataPoints(int numDataPoints) {
	numDataPoints_ = numDataPoints;
}


void Scan::setNumScanOrigins(int numScanOrigins) {
	numScanOrigins_ = numScanOrigins;
	scanOriginNums.resize(numScanOrigins, 0);
	scanOriginParentFileIDs.resize(numScanOrigins, "");
}


Scan::Scan() {
	numDataPoints_ = 0;

	msLevel_ = 0;
	charge_ = -1;
	startMZ_ = -1;
	endMZ_ = -1;
	minObservedMZ_ = -1;
	maxObservedMZ_ = -1;
	basePeakMZ_ = -1;
	basePeakIntensity_ = -1;
	totalIonCurrent_ = -1;
	retentionTimeInSec_ = -1;
	polarity_ = POLARITY_UNDEF;
	analyzer_ = ANALYZER_UNDEF;
	ionization_ = IONIZATION_UNDEF;
	scanType_ = SCAN_UNDEF;
	activation_ = ACTIVATION_UNDEF;

	isCentroided_ = false;
	precursorScanNumber_ = -1;
	precursorScanMSLevel_ = -1;
	precursorCharge_ = -1;
	precursorMZ_ = -1;
	accuratePrecursorMZ_ = false;
	precursorIntensity_ = -1;
	collisionEnergy_ = -1;

	// thermo scans only:
	isThermo_ = false;
	thermoFilterLine_ = "";

	// MassLynx scans only
	isMassLynx_ = false;
	isCalibrated_ = false;

	isMerged_ = false;
	mergedScanNum_ = -1;

	isThresholded_ = false;
	threshold_ = -1;

	// initialize to NULL so that delete can be called alway
	mzArray_ = NULL;
	intensityArray_ = NULL;

	baseline_=NULL;
	noise_=NULL;

	numScanOrigins_ = 0;

	Signal_Th=1e-4;//always using the default value
	RSlower=0.8;
	IsCalibrate=false;
}



Scan::Scan(const Scan& copy) {
	msLevel_ = copy.msLevel_;
	charge_ = copy.charge_;
	startMZ_ = copy.startMZ_;
	endMZ_ = copy.endMZ_;
	minObservedMZ_ = copy.minObservedMZ_;
	maxObservedMZ_ = copy.maxObservedMZ_;
	basePeakMZ_ = copy.basePeakMZ_;
	basePeakIntensity_ = copy.basePeakIntensity_;
	totalIonCurrent_ = copy.totalIonCurrent_;
	retentionTimeInSec_ = copy.retentionTimeInSec_;
	polarity_ = copy.polarity_;
	analyzer_ = copy.analyzer_;
	ionization_ = copy.ionization_;
	scanType_ = copy.scanType_;
	activation_ = copy.activation_;
	isCentroided_ = copy.isCentroided_;
	precursorScanNumber_ = copy.precursorScanNumber_;
	precursorScanMSLevel_ = copy.precursorScanMSLevel_;
	precursorCharge_ = copy.precursorCharge_;
	precursorMZ_ = copy.precursorMZ_;
	accuratePrecursorMZ_ = copy.accuratePrecursorMZ_;
	precursorIntensity_ = copy.precursorIntensity_;
	collisionEnergy_ = copy.collisionEnergy_;
	isThermo_ = copy.isThermo_;
	thermoFilterLine_ = copy.thermoFilterLine_;
	isMassLynx_ = copy.isMassLynx_;
	isCalibrated_ = copy.isCalibrated_;
	isMerged_ = copy.isMerged_;
	mergedScanNum_ = copy.mergedScanNum_;
	isThresholded_ = copy.isThresholded_;
	threshold_ = copy.threshold_;
	

	numScanOrigins_ = copy.numScanOrigins_;

	setNumDataPoints(copy.getNumDataPoints());
	for (int i=0; i<numDataPoints_; i++) {
		mzArray_[i] = copy.mzArray_[i];
		intensityArray_[i] = copy.intensityArray_[i];
		baseline_[i]=copy.baseline_[i];
		noise_[i]=copy.noise_[i];

	}

	//added on 2010.31 _begin
	IsolationWidth_=copy.IsolationWidth_;
	status_par_name.clear();
	size_t i, sn=copy.status_par_name.size();
	for(i=0;i<sn;i++)
	{
		status_par_name.push_back(copy.status_par_name[i]);
	}

	status_par_value.clear();
	sn=copy.status_par_value.size();
	for(i=0;i<sn;i++)
	{
		status_par_value.push_back(copy.status_par_value[i]);
	}

	tailer_par_name.clear();
	sn=copy.tailer_par_name.size();
	for(i=0;i<sn;i++)
	{
		tailer_par_name.push_back(copy.tailer_par_name[i]);
	}

	tailer_par_value.clear();
	sn=copy.tailer_par_value.size();
	for(i=0;i<sn;i++)
	{
		tailer_par_value.push_back(copy.tailer_par_value[i]);
	}

	curScanNum=copy.curScanNum;
	//added on 2010.12.31 section _end
	nativeScanRef_.coordinateType_ = copy.nativeScanRef_.coordinateType_;
	nativeScanRef_.coordinates_.clear();
	for (std::vector<NativeScanRef::CoordinateNameValue>::size_type c = 0;
		c < copy.nativeScanRef_.coordinates_.size();
		c++) {
			nativeScanRef_.coordinates_.push_back(copy.nativeScanRef_.coordinates_[c]);
	}

	//added on 2011.9.5
	IsCalibrate=copy.IsCalibrate;
}



Scan::~Scan() {
	if(numDataPoints_>0) //modified by zhangjy  to make it more safe
	{
		delete [] mzArray_;
		delete [] intensityArray_;

		delete []noise_;
		delete []baseline_;
	}

}


void Scan::centroid(string instrument) {

	// presumed resolution - should be conservative?
	double res400 = 10000.0; // for TOF
	if (instrument == "FT") {
		res400 = 100000.0;
	} else if (instrument == "Orbitrap") {
		res400 = 50000.0;
	}

	// assume sorted
	//  sort(m_peaks.begin(), m_peaks.end(), sortPeaksByMzAsc);

	// determine smallest m/z-interval between peaks
	// this should tell us the frequency with which the
	// mass spectrometer takes readings
	// 
	// this step is helpful because in many profile spectra,
	// the peak is omitted completely if the intensity is below
	// a certain threshold. So the peak list has jumps in m/z values,
	// and neighboring peaks are not necessarily close in m/z. 

	double minInterval = 1000000.0;
	double minIntensity = 1000000.0;

	for (int p = 0; p < numDataPoints_ - 1; p++) {
		//double interval = mzArray_[p + 1] - intensityArray_[p];
        double interval = mzArray_[p + 1] - mzArray_[p];  // fixed typo - DT
		if (interval < 0.0) {
			cerr << "Peak list not sorted by m/z. No centroiding done." << endl;
			return;
		}
		if (minInterval > interval) minInterval = interval;
		if (intensityArray_[p] > 0 /* <- added to ignore zeroes -- DT*/ && minIntensity > intensityArray_[p]) minIntensity = intensityArray_[p];
	}

	vector<Peak>* allPeaks = new vector<Peak>;
	for (int i = 0; i < numDataPoints_ - 1; i++) {
		Peak p;
		p.mz = mzArray_[i];
		p.intensity = (float)intensityArray_[i];
		allPeaks->push_back(p);
		double gap = mzArray_[i + 1] - mzArray_[i];
		double curMz = mzArray_[i];
		int numZeros = 0;
		while (gap > 1.9 * minInterval) {
			if (numZeros < 3 || curMz > mzArray_[i + 1] - 3.1 * minInterval) {
				curMz += minInterval;
			} else {
				curMz = mzArray_[i + 1] - 3.0 * minInterval;
				gap = 4.0 * minInterval;
			}
			Peak pk;
			pk.mz = curMz;
			pk.intensity = 0.0;
			allPeaks->push_back(pk);
			gap -= minInterval;
			numZeros++;
		}

	}


	vector<Peak>* smoothedPeaks = new vector<Peak>;
	for (int i = 0; i < (int)allPeaks->size(); i++) {

		int weight = 6;
		Peak smoothPeak; 
		smoothPeak.mz = (*allPeaks)[i].mz;
		smoothPeak.intensity = 6 * (*allPeaks)[i].intensity;

		if (i >= 2){
			weight += 1;
			smoothPeak.intensity += (*allPeaks)[i - 2].intensity;
		}
		if (i >= 1) {
			weight += 4;
			smoothPeak.intensity += 4 * (*allPeaks)[i - 1].intensity;
		}
		if (i < (int)allPeaks->size() /*numDataPoints_*/ - 1) {   // DT
			weight += 4;
			smoothPeak.intensity += 4 * (*allPeaks)[i + 1].intensity;
		}
		if (i < (int)allPeaks->size() /*numDataPoints_*/ - 2) {   // DT
			weight += 1;
			smoothPeak.intensity += (*allPeaks)[i + 2].intensity;
		}

		smoothPeak.intensity /= (float)weight; 

		smoothedPeaks->push_back(smoothPeak);

	}

	delete (allPeaks);


	//  m_peaks.clear();
	// for (vector<Peak>::iterator sm = smoothedPeaks.begin(); sm != smoothedPeaks.end(); sm++) {
	//   m_peaks.push_back(*sm);
	// }
	// plot(m_pep->getSafeName() + "_smoothed", "");

	int j;
	double maxIntensity;
	int bestPeak;
	bool bLastPos;

	int nextBest;
	double FWHM;

	bLastPos=false;

	vector<Peak> newPeaks;
	//step along each point in spectrum
	for (int i = 0; i < (int)smoothedPeaks->size() - 1; i++){

		//check for change in direction
		if ((*smoothedPeaks)[i].intensity < (*smoothedPeaks)[i + 1].intensity) {

			bLastPos=true;
			continue;

		} else {

			if (bLastPos) {
				bLastPos=false;

				//find max 
				//This is an artifact of using a window of length n (user variable) for identifying
				//a better peak apex on low-res data. Feel free to remove this section if desired.
				//Replace with:
				//  bestPeak=j;
				//  maxIntensity=s[j].intensity;

				maxIntensity = 0;
				for (j = i; j < i + 1; j++) {
					if ((*smoothedPeaks)[j].intensity > maxIntensity) {
						maxIntensity = (*smoothedPeaks)[j].intensity;
						bestPeak = j;
					}
				}

				//Best estimate of Gaussian centroid
				//Get 2nd highest point of peak
				if (bestPeak == smoothedPeaks->size() - 1) {
					nextBest = bestPeak - 1;
				} else if ((*smoothedPeaks)[bestPeak - 1].intensity > (*smoothedPeaks)[bestPeak + 1].intensity) {
					nextBest = bestPeak - 1;
				} else {
					nextBest = bestPeak + 1;
				}

				//Get FWHM of Orbitrap
				//This is the function you must change for each type of instrument.
				//For example, the FT would be:
				//  FWHM = s[bestPeak].mz * s[bestPeak].mz / (400*res400);
				// Orbitrap
				// FWHM = (*smoothedPeaks)[bestPeak].mz * sqrt((*smoothedPeaks)[bestPeak].mz) / (20 * res400);

				if (instrument == "FT") {
					FWHM = (*smoothedPeaks)[bestPeak].mz * (*smoothedPeaks)[bestPeak].mz / (400 * res400);
				} else if (instrument == "Orbitrap") {
					FWHM = (*smoothedPeaks)[bestPeak].mz * sqrt((*smoothedPeaks)[bestPeak].mz) / (20 * res400);
				} else {
					// for TOF
					FWHM = (*smoothedPeaks)[bestPeak].mz / res400;
				}

				Peak centroid;
				//Calc centroid MZ (in three lines for easy reading)
				centroid.mz = pow(FWHM , 2) * log((*smoothedPeaks)[bestPeak].intensity / (*smoothedPeaks)[nextBest].intensity);
				centroid.mz /= 8 * log(2.0) * ((*smoothedPeaks)[bestPeak].mz - (*smoothedPeaks)[nextBest].mz);
                if ( fabs(centroid.mz) < fabs( ((*smoothedPeaks)[bestPeak].mz - (*smoothedPeaks)[nextBest].mz) / 2 ) ) // sanity check - DT
				    centroid.mz += ((*smoothedPeaks)[bestPeak].mz + (*smoothedPeaks)[nextBest].mz) / 2;
                else
                    centroid.mz = (*smoothedPeaks)[bestPeak].mz; // fail-safe - DT


				//Calc centroid intensity
				// double exponent = pow(((*smoothedPeaks)[bestPeak].mz - centroid.mz) / FWHM, 2) * (4 * log(2.0));
				//
				// if (exponent > 1.0) {
				//  centroid.intensity = (*smoothedPeaks)[bestPeak].intensity;
				// } else {
				//  centroid.intensity = (*smoothedPeaks)[bestPeak].intensity / exp(-exponent);
				// }

				centroid.intensity = (*smoothedPeaks)[bestPeak].intensity;

				//Hack until I put in mass ranges
				//Another fail-safe can be made for inappropriate intensities
				if (centroid.mz < 0 || centroid.mz > 2000 || centroid.intensity < 0.99 * minIntensity) {
					//do nothing if invalid mz
				} else {
					newPeaks.push_back(centroid);
				}

			}

		}
	}

	delete (smoothedPeaks);

	vector<Peak> finalPeakList;
	finalPeakList.clear();
	// m_origMaxIntensity = 0.0;
	// m_totalIonCurrent = 0.0;
	// m_isAnnotated = false;

	// reset Scan data arrays
	setNumDataPoints((int)newPeaks.size());

	// copy centroided results to Scan data arrays
	long ci=0;
	for (vector<Peak>::iterator j = newPeaks.begin(); j != newPeaks.end(); j++) {
		mzArray_[ci] = (*j).mz;
		intensityArray_[ci] = (*j).intensity;
        ci ++; // added -- DT
		//    if (m_origMaxIntensity < j->intensity) m_origMaxIntensity = j->intensity;
		//    m_totalIonCurrent += j->intensity;
	}

	isCentroided_ = true;

	// TODO: reset/recalc other scan values? 

	// m_isScaled = false;

	// if (m_bins) delete (m_bins);
	// m_bins = false;

	// delete (m_intensityRanked);
	// m_intensityRanked = NULL;

	//  plot(m_pep->getSafeName() + "_centroided", "");
}


// if not discard, rewrite as zero
void Scan::threshold(double inclusiveCutoff, bool discard) {

	vector<Peak> newPeakList;
	newPeakList.clear();

	int i;
	int orig = numDataPoints_;

	/*for (i=0; i<numDataPoints_; i++) {
		cout << mzArray_[i] << "\t" << intensityArray_[i] << endl;
	}
	cout << endl << endl;*/


	for (i=0; i<numDataPoints_; i++) {
		double curIntensity = intensityArray_[i];
		double curMZ = mzArray_[i];
		if (curIntensity >= inclusiveCutoff) {
			// save the value
			Peak p;
			p.intensity = curIntensity;
			p.mz = curMZ;
			newPeakList.push_back(p);
		}
		else {
			// this failed the cutoff.
			// do we discard it, or zero it out?
			if (!discard) {
					// insert a zero'd intensity here
					Peak p;
					p.intensity = 0;
					p.mz = curMZ;
					newPeakList.push_back(p);
			}
			/*else {
				cout << "discarded " << curMZ << "\t" << curIntensity << endl;
			}*/
		}
	}

	// rewrite our actual internal data:

	// reset Scan data arrays:
	setNumDataPoints((int)newPeakList.size());
	//if (numDataPoints_ > 0){
	//	//cout << "scan " <<  
	//	cout << "threshold: cutoff is " << inclusiveCutoff
	//	<< " discard is " << discard << endl;
	//	cout << "original #data points: " << numDataPoints_ << endl;
	//	cout << "new #data points: " << numDataPoints_ << endl;
	//}
	//cout << endl << endl;

	// copy thresholded results to Scan data arrays
	long ci=0;
	for (vector<Peak>::iterator j = newPeakList.begin(); j != newPeakList.end(); j++) {
		mzArray_[ci] = (*j).mz;
		intensityArray_[ci] = (*j).intensity;
		ci++;
	}

	isThresholded_ = true;

	// TODO: recalulate other Scan values (basepeak, range, etc)?
}

double Scan::GetSampleDm(std::string instrument,double mass)
{
	if(instrument=="FT") return 6.06e-9*mass*mass-7.063e-11*mass+4.085e-8;//for FT
	else if(instrument=="Orbitrap") return 1.12e-009*mass*mass+2.204e-006*mass-0.0003366;//for Orbitrap
	else return 0;
}

double Scan::GetFitInt(double cen_mz,string instrument)
{
	if(numDataPoints_<=0) return 0;
	if(isCentroided_) return 0;
	double dm;
	dm=1.5*GetSampleDm(instrument,cen_mz);
	int middle=TFindL(0,numDataPoints_,cen_mz);
	if(middle==-1) return 0;
	int begin=middle-1;
	double mzold=mzArray_[middle];
	while(begin>=0)
	{
		if(intensityArray_[begin]<1e-2) break;
		if(mzArray_[begin]+dm<mzold) break;
		mzold=begin;
		begin--;
	}
	if(begin<0) begin=0;
	int end=middle+1;
	mzold=mzArray_[middle];
	while(end<numDataPoints_)
	{
		if(intensityArray_[end]<1e-2) break;
		if(mzArray_[end]>mzold+dm) break;
		mzold=end;
		end++;
	}
	if(end>=numDataPoints_) end=numDataPoints_-1;	
	int i;
	vector<double> pkmass;
	vector<double> pkintensity;
	vector<double> beta;
	for(i=begin;i<end;i++)
	{	
		pkmass.push_back(mzArray_[i]);
		pkintensity.push_back(intensityArray_[i]);
	}
	double RSquare=Gauss_Fit(pkmass,pkintensity,beta);
	int n=beta.size()/4;	
	double DM=1.0;
	int f_idx=-1;
	for(i=0;i<n;i++)
	{
		double dm=fabs(cen_mz-beta[4*i+1]);
		if(dm<DM) 
		{
			DM=dm;
			f_idx=i;
		}
	}
	if(f_idx==-1) return 0;
	else return  beta[4*f_idx]*beta[4*f_idx+2]*SQRTPI;
}

double Scan::GetFitInt(double cen_mz)
{
	if(numDataPoints_<=0) return 0;
	if(isCentroided_) return 0;
	string instrument="FT";
	double dm;
	dm=1.5*GetSampleDm(instrument,cen_mz);
	int middle=TFindL(0,numDataPoints_,cen_mz);
	if(middle==-1) return 0;
	int begin=middle-1;
	double mzold=mzArray_[middle];
	while(begin>=0)
	{
		if(intensityArray_[begin]<1e-2) break;
		if(mzArray_[begin]+dm<mzold) break;
		mzold=begin;
		begin--;
	}
	if(begin<0) begin=0;
	int end=middle+1;
	mzold=mzArray_[middle];
	while(end<numDataPoints_)
	{
		if(intensityArray_[end]<1e-2) break;
		if(mzArray_[end]>mzold+dm) break;
		mzold=end;
		end++;
	}
	if(end>=numDataPoints_) end=numDataPoints_-1;	
	int i;
	vector<double> pkmass;
	vector<double> pkintensity;
	vector<double> beta;
	for(i=begin;i<=end;i++)
	{	
		pkmass.push_back(mzArray_[i]);
		pkintensity.push_back(intensityArray_[i]);
	}
	double RSquare=Gauss_Fit(pkmass,pkintensity,beta);
	int n=beta.size()/4;	
	double DM=1.0;
	int f_idx=-1;
	for(i=0;i<n;i++)
	{
		double dm=fabs(cen_mz-beta[4*i+1]);
		if(dm<DM) 
		{
			DM=dm;
			f_idx=i;
		}
	}
	if(f_idx==-1) return 0;
	else return  beta[4*f_idx]*beta[4*f_idx+2]*SQRTPI;
}

int Scan::Centriod_gf(std::string instrument)
{	
	if(numDataPoints_<=0) return 0;
	if(isCentroided_) return 0;

	double dm;
	int begin=0;	
	vector<double> c_intarray;
	vector<double> c_mzarray;
	int i;
	while(begin<numDataPoints_)
	{				
		int end=begin+1;
		int idx_old=begin;
		dm=1.5*GetSampleDm(instrument,mzArray_[begin]);
		while(end<numDataPoints_)
		{
			if(intensityArray_[end]<Signal_Th) break;
			if(mzArray_[end]>mzArray_[idx_old]+dm) break;			
			idx_old=end;
			end++;
		}
		vector<double> pkmass;
		vector<double> pkintensity;
		vector<double> beta;
		for(i=begin;i<end;i++)
		{	
			pkmass.push_back(mzArray_[i]);
			pkintensity.push_back(intensityArray_[i]);
		}
		double RSquare=Gauss_Fit(pkmass,pkintensity,beta);
		if(RSquare>RSlower)
		{
			int n=beta.size()/4;
			int mu_idx;
			int abu_idx;
			for(i=0;i<n;i++)
			{
				abu_idx=4*i;
				mu_idx=abu_idx+1;
				if(beta[abu_idx]<=0||beta[mu_idx]<=0) continue;
				c_mzarray.push_back(beta[mu_idx]);
				c_intarray.push_back(beta[abu_idx]);
			}
		}
		begin=end;
	}

	delete []mzArray_;
	delete []intensityArray_;
	numDataPoints_=c_intarray.size();
	if(numDataPoints_<=0) return 0;
	mzArray_=new double[numDataPoints_];
	intensityArray_=new double[numDataPoints_];
	for(i=0;i<numDataPoints_;i++)
	{
		intensityArray_[i]=c_intarray[i];
		mzArray_[i]=c_mzarray[i];
	}	
	isCentroided_ = true;
	return numDataPoints_;
}

bool Scan::OutPutData(FILE *fp)
{
	if(fp==NULL) return false;
	//fprintf(fp,"%d.%d\n",curScanNum,numDataPoints_);
	for(int i=0;i<numDataPoints_;i++)
	{
		fprintf(fp,"%lf\t%lf\n",mzArray_[i],intensityArray_[i]);
	}
	return true;
}

int Scan::TFindL(int b,int e,double fmz)
{
	if(e>=numDataPoints_) e=numDataPoints_-1;
	if(b>=e) return -1;
	if(mzArray_[b]>fmz+0.01) return -1;
	if(mzArray_[e]<fmz-0.01) return -1;
	if(e-b<=1) return b;	
	int idx=(b+e)/2;
	if(mzArray_[idx]>=fmz)
	{
		return TFindL(b,idx,fmz);
	}
	else return TFindL(idx,e,fmz);
}

int Scan::TFindH(int b,int e,double fmz)
{
	if(e>=numDataPoints_) e=numDataPoints_-1;
	if(b>=e) return -1;
	if(mzArray_[b]>fmz+0.01) return -1;
	if(mzArray_[e]<fmz-0.01) return -1;
	if(e-b<=1) return e;	
	int idx=(b+e)/2;
	if(mzArray_[idx]>=fmz)
	{
		return TFindL(b,idx,fmz);
	}
	else return TFindL(idx,e,fmz);
}

int Scan::FindByPmz(double pmz,double dm)
{
	if(numDataPoints_<=0) return -1;
	int begin=TFindL(0,numDataPoints_-1,pmz-dm);
	if(begin==-1) return -1;
	int end=TFindH(0,numDataPoints_-1,pmz+dm);
	if(end==-1) return -1;
	double MaxInt=0;
	int final=-1;
	while(begin<numDataPoints_&&mzArray_[begin]<pmz-dm) begin++;
	while(end>=0&&mzArray_[end]>pmz+dm) end--;	
	for(int i=begin;i<=end;i++)
	{	
		if(MaxInt<intensityArray_[i]) 
		{
			MaxInt=intensityArray_[i];
			final=i;
		}
	}	
	return final;
}


bool Scan::GetPosChForMS2(vector<int> CH)
{
	if(msLevel_!=2) return false;
	if(numDataPoints_<2) return false;
	if(precursorCharge_>0) return false;
	double pi[7];
	int i;
	for(i=0;i<7;i++)pi[i]=0;
	for(i=1;i<numDataPoints_;i++)
	{	
		pi[5]+=intensityArray_[i];
		if(mzArray_[i]<0.5f*precursorMZ_) pi[0]+=intensityArray_[i];
		else if(mzArray_[i]<precursorMZ_) pi[1]+=intensityArray_[i];
		else if(mzArray_[i]<1.5f*precursorMZ_) pi[2]+=intensityArray_[i];
		else if(mzArray_[i]<2.0f*precursorMZ_) pi[3]+=intensityArray_[i];
		else pi[4]+=intensityArray_[i];
	}		
	pi[6]=endMZ_-startMZ_;
	//double B[8]={7543.50211,7545.88562,7547.90369,7546.62426,7609.729,.000839975,-0.000311281,-7544.3281};
	double BLTQ[8]={-34.1998475403346,-31.2367136981516,-30.4033514691468,-30.5198459718126,0,0.000158113380637447,0.118210010377885,32.8992182840002};
	double chidx=0;
	if(pi[5]<=1e-4) return false;
	for(i=0;i<5;i++) chidx+=BLTQ[i]*pi[i]/pi[5];
	chidx+=BLTQ[5]*(pi[5])+BLTQ[7]+BLTQ[6]*pi[6];
	if(chidx<1.5) CH.push_back(1);
	else if(chidx<2.0) CH.push_back(2);
	else if(chidx<3.0) 
	{
		CH.push_back(2);
		CH.push_back(3);
	}
	else CH.push_back(3);
	return true;	
}
