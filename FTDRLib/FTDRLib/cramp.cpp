/***************************************************************************
cramp.cpp

/***************************************************************************
cramp.hpp

  A C++ wrapper for the RAMP code.

  Use this library to parse an mzXML file in a non-sequential way, by
  taking advantage of the index element.  

  (C) 2004 by Brian Pratt, Insilicos LLC 

  Based on mzXML2Other, which has this copyright:
	 -------------------
	 begin					 : Wed Apr 2
	 copyright				 : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
	 email					 : ppatrick@student.ethz.ch
 ***************************************************************************/

/***************************************************************************
*																								  *
*	 This program is free software; you can redistribute it and/or modify  *
*	 it under the terms of the GNU Library or "Lesser" General Public 	  *
*	 License (LGPL) as published by the Free Software Foundation;			  *
*	 either version 2 of the License, or (at your option) any later		  *
*	 version.																				  *
*																								  *
***************************************************************************/



#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "stdio.h"
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include "sys/errno.h"
#endif
#include "cramp.hpp"

#include "math.h"



/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file.
 *
 * @param fileName: Name of the msxml file
 * @param startSCan: Number of the scan we want to read from
 * @param what: -HEADER will return num, msLevel and retentionTime
 *              -SCAN will return only the peaks
 *              -ALL will return everything found in scan, precursorMz and peaks
 *
 * @return pData is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!!
 */



cRamp::cRamp( const char* fileName,bool declaredScansOnly ) : 
  m_filename(fileName), m_declaredScansOnly(declaredScansOnly), m_runInfo()
{
   m_handle = rampOpenFile(fileName);
   m_scanOffsets = NULL;
   m_runInfo = NULL;
   m_lastScan = 0;   
   if (!OK()) {
     // HENRY -- I would prefer this to be silent, and let the caller deals with it
     // cout << "Error: Could not open file " << fileName << ": " << strerror(errno) << endl;
     // END HENRY
   } else {

      m_runInfo = getRunInfo();

      // HENRY -- always read index to set scan count, since scan count
      // declared at the top of the mzXML file is unreliable now that
      // there are missing scans.
      // This will also set the structs m_scanOffsets, and the value m_lastScan
      
      //      if (m_runInfo->m_data.scanCount < 0) { // undeclared scan count
         // this will provoke reading of index, which sets scan count
         rampScanInfo* tmp = getScanHeaderInfo ( 1 );
         free(tmp);
	 // }
      // END HENRY
   }
}

cRamp::~cRamp() {
   rampCloseFile(m_handle);
   
   // NB: these pointers may be null on file open failure, 
   // but free/delete of NULL is OK per C++ standard 
   free(m_scanOffsets);
   delete m_runInfo; // was free() - but allocated with new 
}

//
// here are the private guts
//
rampInfo* cRamp::do_ramp( ramp_fileoffset_t arg , eWhatToRead	what )
{
   
   switch( what ) {
   case RAMP_RUNINFO:
   case RAMP_HEADER:
   case RAMP_PEAKS:
   case RAMP_INSTRUMENT:
      break; // OK
   default:
	  std::cerr << "unknown read type!\n";
      return NULL;
      break;
   }	
   
   rampInfo* returnPtr=NULL;
   
   if ((RAMP_RUNINFO != what) && (RAMP_INSTRUMENT != what) && !m_scanOffsets) {
      int iLastScan = 0; 
     // we need the index to get anything besides the header
      ramp_fileoffset_t indexOffset = getIndexOffset(m_handle);
      m_scanOffsets = readIndex(m_handle, indexOffset, &iLastScan);
      if (iLastScan >= m_runInfo->m_data.scanCount) {
		 if (!m_declaredScansOnly) {
         m_runInfo->m_data.scanCount = iLastScan;
		 } else { // get rid of all the fake entries created
			 for (int n=1;n<=iLastScan;n++) { // ramp is 1 based
				 if (m_scanOffsets[n]==-1) {
					// find a run of fakes
				    int m;
					for (m=n+1;(m<=iLastScan)&&(m_scanOffsets[m]==-1);m++);
					if (m<=iLastScan) {
						memmove(m_scanOffsets+n,m_scanOffsets+m,
						  sizeof(ramp_fileoffset_t)*((iLastScan-m)+1));
      }
					iLastScan-=(m-n);
				 }
			 }
		 }
      }
      // HENRY - store last scan explicitly.
      m_lastScan = iLastScan;
      // END HENRY
   }

   
   // HENRY -- arg is out of bounds. instead of creating havoc in RAMP, let's just kill it here.
   if (RAMP_RUNINFO != what && (RAMP_INSTRUMENT != what) && (arg > m_runInfo->m_data.scanCount || arg < 1)) {
     return (NULL);
   }
     
   if (m_scanOffsets || (RAMP_RUNINFO == what) || (RAMP_INSTRUMENT == what)) {
      ramp_fileoffset_t scanOffset=-1;
      if (RAMP_RUNINFO == what || RAMP_INSTRUMENT == what) {
         scanOffset = 0; // read from head of file
      } else {
         scanOffset = m_scanOffsets[arg]; // ramp is one-based
      }
      
      if (scanOffset >= 0) {
         
         // -----------------------------------------------------------------------
         // And now we can parse the info we were looking for
         // -----------------------------------------------------------------------
         
         
         // Ok now we have to copy everything in our structure
         switch( what )
         {
         case RAMP_RUNINFO:
            returnPtr = new rampRunInfo( m_handle );
            break;
         case RAMP_HEADER:
            returnPtr = new rampScanInfo( m_handle, scanOffset, (int)arg );
            if (returnPtr) {
#ifdef HAVE_PWIZ_MZML_LIB
			   if (!m_handle->mzML) // rampadapter already set this for us
#endif
              ((rampScanInfo *)returnPtr)->m_data.filePosition = scanOffset; // for future reference
            
              // HENRY -- error checking here
              if (((rampScanInfo*)returnPtr)->m_data.acquisitionNum < 0) {
                // something failed in RAMP, possibly because it's a missing scan
                delete ((rampScanInfo*)returnPtr);
                returnPtr = NULL;
              }
            }
            break;           
         case RAMP_PEAKS:
            returnPtr = new rampPeakList( m_handle, scanOffset);
            
            // HENRY -- error checking here
            if (returnPtr && ((rampPeakList*)returnPtr)->getPeakCount() <= 0) {
              // something failed in RAMP, possibly because it's a missing scan
              delete ((rampPeakList*)returnPtr);
              returnPtr = NULL;
            }
            break;
            
         // HENRY -- add the instrument info reading functionality (present in RAMP, but not provided in cRAMP before)
         case RAMP_INSTRUMENT:
            returnPtr = new rampInstrumentInfo(m_handle);
            if (((rampInstrumentInfo*)returnPtr)->m_instrumentStructPtr == NULL) {
              delete ((rampInstrumentInfo*)returnPtr);
              returnPtr = NULL;
            }
            break;
         }
         
      }
   }
   
   
   
   return returnPtr;
}


/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file to obtain minimal info on the msRun contained in the file.  

 *
 * @return rapRunInfo* is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!!
 */

rampRunInfo* cRamp::getRunInfo (  ) {
   rampRunInfo* result;
   if (m_runInfo) { // did we derive this already?
      result = new rampRunInfo(*m_runInfo);
   } else {
      result = (rampRunInfo*) do_ramp(0, RAMP_RUNINFO);
   }
   return result;
}

/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file to obtain minimal header info for a numbered scan (thus minimizing parse time).  
 *
 * @param fileName: Name of the msxml file
 * @param startSCan: Number of the scan we want to read from
 * @return rapHeaderInfo* is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!! returns just the minimal header info num, msLevel and retentionTime
 */

rampScanInfo* cRamp::getScanHeaderInfo ( int whichScan  ) {
   return (rampScanInfo*) do_ramp((ramp_fileoffset_t)whichScan, RAMP_HEADER);
}


/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file to obtain peak info for a numbered scan.
 *
 * @param fileName: Name of the msxml file
 * @param startSCan: Number of the scan we want to read from
 * @return rapPeakList* is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!! returns everything found in scan, precursorMz and peaks
 */

rampPeakList* cRamp::getPeakList ( int whichScan ) {
   return (rampPeakList*) do_ramp((ramp_fileoffset_t)whichScan, RAMP_PEAKS);
}

// HENRY - provides instrument info getting method
rampInstrumentInfo* cRamp::getInstrumentInfo () {
  return (rampInstrumentInfo*) do_ramp(0, RAMP_INSTRUMENT);
}
// END HENRY


// HENRY - sequential access parser that skips over missing scans. This version only reads scan header.
bool cRampIterator::nextScan(rampScanInfo** scanInfo) {
	while (++m_currentScan <= m_cramp.getLastScan() && m_cramp.getScanOffset(m_currentScan) <= 0);
  if (m_currentScan > m_cramp.getLastScan()) {
    return (false);
  }
  
  *scanInfo = (rampScanInfo*)m_cramp.do_ramp((ramp_fileoffset_t)(m_currentScan), RAMP_HEADER);
  return (true);
}
// END HENRY

// HENRY - sequential access parser that skips over missing scans. This version reads both scan header and peak list.
bool cRampIterator::nextScan(rampScanInfo** scanInfo, rampPeakList** peakList) {
  while (++m_currentScan <= m_cramp.getLastScan() && m_cramp.getScanOffset(m_currentScan) <= 0);
  if (m_currentScan > m_cramp.getLastScan()) {
    return (false);
  }
  
  *scanInfo = (rampScanInfo*)m_cramp.do_ramp((ramp_fileoffset_t)(m_currentScan), RAMP_HEADER);
  *peakList = (rampPeakList*)m_cramp.do_ramp((ramp_fileoffset_t)(m_currentScan), RAMP_PEAKS);

  return (true);

}
// END HENRY

// HENRY - resets the sequential access parser to the first scan.
void cRampIterator::reset() {
  m_currentScan = 1;
}
// END HENRY

/**
 * populate from a file handle
 **/
rampPeakList::rampPeakList(RAMPFILE *handle, ramp_fileoffset_t index) {
   init();
   m_peaksCount = readPeaksCount(handle,index);
   m_pPeaks = (rampPeakInfoStruct *)readPeaks(handle,index);
}

/**
 * populate from a file handle
 **/
rampScanInfo::rampScanInfo(RAMPFILE *handle, ramp_fileoffset_t index, int seqNum) {
   init();
   readHeader(handle,index,&m_data);
   m_data.seqNum = seqNum;
}

/**
 * populate from a file handle
 **/
rampRunInfo::rampRunInfo(RAMPFILE *handle) {
   init();
   readMSRun(handle,&m_data);
}

// HENRY - provides instrument info reading functionality
rampInstrumentInfo::rampInstrumentInfo(RAMPFILE *handle) {
  init();
  m_instrumentStructPtr = getInstrumentStruct(handle);
}
//add by zhangjy 2008.5.30
bool cRamp::GetTIC(int ScanNum,double *TIC, double *RT)
{
	*TIC=0;
	*RT=0;
	int Last=getLastScan();
	int First=1;
	rampScanInfo *pt;
	if(ScanNum<First||ScanNum>Last) return false;
	pt=getScanHeaderInfo(ScanNum);
	int tscan=ScanNum;
	while(pt->m_data.msLevel>1) 
	{
		tscan--;
		if(tscan<First) return false;
		pt=getScanHeaderInfo(tscan);
	}
	*RT=pt->m_data.retentionTime;
	*TIC=pt->m_data.totIonCurrent;	
	delete pt;
	return true;
}
//added by zhangjy on 2010.1.31 for the global calibration
//call this function before any operation after the object creation
void cRamp::Initial()
{
	charge=0;
	bMass=0;
	pkl=NULL;
	sFName[0]='\0';
	Last=getLastScan();
	First=1;
	MET=20;
}

double cRamp::match(mPeak ISOP[MAX_ISO],double IsoDisT[MAX_ISO])
{
	double sumx,sumy,sumxy;
	sumx=0;
	sumy=0;
	sumxy=0;
	for(int i=0;i<MAX_ISO;i++)
	{
		sumx+=ISOP[i].dIntensity*ISOP[i].dIntensity;
		sumy+=IsoDisT[i]*IsoDisT[i];
		sumxy+=ISOP[i].dIntensity*IsoDisT[i];
	}	
	sumx*=sumy;
	sumx=sqrt(sumx);
	if(sumx<=0) return 0;	
	return sumxy/sumx;
}

bool cRamp::PreFindMInt(IsoCluster &it)
{
	if(charge<=0) return false;	
	double dm=0.05;//charge*bMass*MET/1e6;//20ppm
	double step=1.007825/charge;
	double subLP[MAX_ISO+3],subUP[MAX_ISO+3];
	mPeak tmpISO[MAX_ISO+3];
	long i,k;
	for(i=0;i<MAX_ISO+3;i++)
	{
		tmpISO[i].dIntensity=0;
		double Cen=bMass+(i-3)*step;
		subLP[i]=Cen-dm;
		subUP[i]=Cen+dm;
	}
	
	double basePeak=0;
	int PeakCount=pkl->getPeakCount();
	rampPeakInfoStruct *pt;
	for(i=0;i<PeakCount;i++) 
	{
		pt=pkl->getPeak(i);
		if(basePeak<pt->intensity) basePeak=pt->intensity;
	}
	int Begin=TFindL(pkl,0,PeakCount,subLP[0]);
	int End=TFindL(pkl,0,PeakCount,subUP[MAX_ISO-1]);
	if(End<PeakCount-1) End++;
	for(i=Begin;i<End;i++)
	{	
		pt=pkl->getPeak(i);
		for(k=0;k<MAX_ISO+3;k++)
		{
			if(pt->mz>subLP[k]&&pt->mz<subUP[k])
			{
				if(tmpISO[k].dIntensity<pt->intensity) 
				{
					tmpISO[k].dMass=pt->mz;	
					tmpISO[k].dIntensity=pt->intensity;
				}
			}
		}			
	}

	int basei=3;
	double IsoDisT[MAX_ISO];
	CalDefIsoDis(bMass*charge,IsoDisT);
	double gd_max=match(tmpISO+3,IsoDisT);;
	double gd;
	for(i=2;i>=0;i--)
	{
		gd=match(tmpISO+i,IsoDisT);
		if(gd>gd_max) 
		{
			basei=i;
			gd_max=gd;
		}
	}
	if(gd_max<GD_CUT) return false;
	for(i=0;i<MAX_ISO;i++)
	{
		it.ISOP[i]=tmpISO[i+basei];

	}
	if(basePeak>1e-4)
	{
		for(i=0;i<MAX_ISO;i++)
		{
			it.ISOP[i].rInt=it.ISOP[i].dIntensity/basePeak;
		}
	}	
	if(it.ISOP[0].rInt<RINT_CUT) return false;
	it.charge=charge;
	bMass=it.ISOP[0].dMass;
	if(bMass<MASS_LIM) return false;
	return true;
}

bool cRamp::FindMInt(IsoCluster &it)
{
	if(charge<=0) return false;	
	double dm=charge*bMass*MET/1e6;//20ppm
	double step=1.007825/charge;
	double subLP[MAX_ISO],subUP[MAX_ISO];
	long i,k;
	for(i=0;i<MAX_ISO;i++)
	{
		it.ISOP[i].dIntensity=0;
		double Cen=bMass+i*step;
		subLP[i]=Cen-dm;
		subUP[i]=Cen+dm;
	}
	
	double basePeak=0;
	int PeakCount=pkl->getPeakCount();
	rampPeakInfoStruct *pt;
	for(i=0;i<PeakCount;i++) 
	{
		pt=pkl->getPeak(i);
		if(basePeak<pt->intensity) basePeak=pt->intensity;
	}
	int Begin=TFindL(pkl,0,PeakCount,subLP[0]);
	int End=TFindL(pkl,0,PeakCount,subUP[MAX_ISO-1]);
	if(End<PeakCount-1) End++;
	for(i=Begin;i<End;i++)
	{	
		pt=pkl->getPeak(i);
		for(k=0;k<MAX_ISO;k++)
		{
			if(pt->mz>subLP[k]&&pt->mz<subUP[k])
			{
				if(it.ISOP[k].dIntensity<pt->intensity) 
				{
					it.ISOP[k].dMass=pt->mz;	
					it.ISOP[k].dIntensity=pt->intensity;
				}
			}
		}			
	}
	if(basePeak>1e-4)
	{
		for(i=0;i<MAX_ISO;i++)
		{
			it.ISOP[i].rInt=it.ISOP[i].dIntensity/basePeak;
		}
	}
	if(it.ISOP[0].rInt<RINT_CUT) return false;
	double IsoDisT[MAX_ISO];
	CalDefIsoDis(bMass*charge,IsoDisT);
	if(it.Match(IsoDisT)<GD_CUT) return false;
	it.charge=charge;
	return true;
}

int cRamp::TFindL(rampPeakList *pkl,int b,int e,double fmz)
{
	if(e-b<=1) return b;	
	int idx=(b+e)/2;
	rampPeakInfoStruct *pt;
	pt=pkl->getPeak(idx);
	if(pt->mz>=fmz)
	{
		return TFindL(pkl,b,idx,fmz);
	}
	else return TFindL(pkl,idx,e,fmz);
}

bool cRamp::RefinePISO(IsoCluster &it)
{
	if(FindMInt(it)) return true;
	return false;
}

bool cRamp::PreRefinePISO(IsoCluster &it)
{
	for(charge=MAX_CH;charge>=1;charge--)
	{			
		if(PreFindMInt(it)) return true;
	}
	return false;
}

void cRamp::GetIsoCData(IsoCluster *it,vector<double> *dMass,double C[5])
{
	for(int i=0;i<MAX_ISO;i++)
	{
		if(it->ISOP[i].rInt<RINT_CUT) break;
		double mze=it->ISOP[i].dMass/MZE_SCALE;
		double tic=it->pTIC/TIC_SCALE;
		double rt=it->RT;
		double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*rt-ISO_DIFF[i]/charge;
		dMass->push_back(cmass);
	}
}

void cRamp::GetXICData(vector<IsoCluster> *XICList,vector<double> *dMass,double C[5])
{
	dMass->clear();
	size_t i,s=XICList->size();	
	for(i=0;i<s;i++)
	{
		GetIsoCData(&(XICList->at(i)),dMass,C);
	}
}

void cRamp::CalMST(vector<double> *dMass,double *m,double *s)
{
	*m=0;
	*s=0;
	size_t i,size=dMass->size();
	if(size<=0) return;
	for(i=0;i<size;i++)
	{
		*m+=dMass->at(i);
	}
	*m/=size;

	for(i=0;i<size;i++)
	{
		double t=dMass->at(i)-*m;
		*s+=t*t;
	}
	*s=sqrt(*s/size);
}

double cRamp::GetMSTIC()
{
	double pTIC=0;
	int PeakCount=pkl->getPeakCount();
	for(int i=0;i<PeakCount;i++) pTIC+=pkl->getPeak(i)->intensity;
	return pTIC;
}


double cRamp::GetGCError(double C[5], long ScanNum,double pMZ,int ch)
{
	rampScanInfo *pt;
	if(ScanNum<First||ScanNum>Last) return 0;
	pt=getScanHeaderInfo(ScanNum);
	while(pt->m_data.msLevel!=2)  
	{
		delete pt;
		return 0;
	}
	delete pt;
	//refine the parent mz value in the MS scan
	bMass=pMZ;
	charge=ch;
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>First)
	{
		tempSNo--;
		pt=getScanHeaderInfo(tempSNo);
		if(pt->m_data.msLevel==1)  break;
	}
	if(tempSNo<First) return 0;
	pkl=getPeakList(tempSNo);
	//read the MS scan	
	int PeakCount=pkl->getPeakCount();
	if(PeakCount<=0) return 0;
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;
	if(!RefinePISO(it))
	{
		delete pkl;
		return 0;
	}
	it.pTIC=GetMSTIC();
	it.RT=pt->getRetentionTimeSeconds()/60;
	delete pt;
	delete pkl;
	XICList_F.push_back(it);
	int Interupt_times=0;
	while(tempSNo>First)
	{
		tempSNo--;
		pt=getScanHeaderInfo(tempSNo);
		if(pt->m_data.msLevel!=1)  continue;
		pkl=getPeakList(tempSNo);		
		if(pkl->getPeakCount()<=0) 
		{
			delete pt;
			continue;	
		}
		if(RefinePISO(it)) 
		{
			it.pTIC=GetMSTIC();
			it.RT=pt->getRetentionTimeSeconds()/60;
			XICList_F.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;
		delete pkl;
		delete pt;
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
	Interupt_times=0;
	while(tempSNo<Last)
	{
		tempSNo++;
		pt=getScanHeaderInfo(tempSNo);
		if(pt->m_data.msLevel!=1)  continue;
		pkl=getPeakList(tempSNo);
		if(pkl->getPeakCount()<=0) 
		{
			delete pt;
			continue;	
		}
		if(RefinePISO(it)) 
		{
			it.pTIC=GetMSTIC();
			it.RT=pt->getRetentionTimeSeconds()/60;
			XICList_F.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;
		delete pkl;
		delete pt;
		if(Interupt_times>=IT_LIM) break;
	}

	vector<double> dMass;
	GetXICData(&XICList,&dMass,C);
	double m,s;
	CalMST(&dMass,&m,&s);
	return (m-pMZ)*1e6/pMZ;
}

bool cRamp::XICCalibrate(double C[5],FILE *MGFfp,long ScanNum)
{
	rampScanInfo *pt;
	if(ScanNum<First||ScanNum>Last) return 0;
	pt=getScanHeaderInfo(ScanNum);
	while(pt->m_data.msLevel!=2)  
	{
		delete pt;
		return 0;
	}
	delete pt;
	//refine the parent mz value in the MS scan
	
	bMass=pt->m_data.precursorMZ;
	charge=pt->m_data.precursorCharge;
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>First)
	{
		tempSNo--;
		pt=getScanHeaderInfo(tempSNo);
		if(pt->m_data.msLevel==1)  break;
	}
	if(tempSNo<First) return false;
	pkl=getPeakList(tempSNo);
	//read the MS scan	
	if(pkl->getPeakCount()<=0) return false;
	//refine the pmz
	IsoCluster it;
	vector<IsoCluster> XICList_F;
	if(!RefinePISO(it))
	{
		delete pkl;
		return false;
	}
	it.pTIC=GetMSTIC();
	it.RT=pt->getRetentionTimeSeconds()/60;
	delete pt;
	delete pkl;
	XICList_F.push_back(it);
	int Interupt_times=0;
	while(tempSNo>First)
	{
		tempSNo--;
		pt=getScanHeaderInfo(tempSNo);
		if(pt->m_data.msLevel!=1)  continue;
		pkl=getPeakList(tempSNo);
		if(pkl->getPeakCount()<=0) 
		{
			delete pt;
			continue;	
		}
		if(RefinePISO(it)) 
		{
			it.pTIC=GetMSTIC();
			it.RT=pt->getRetentionTimeSeconds()/60;
			XICList_F.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;
		delete pkl;
		delete pt;
		if(Interupt_times>=IT_LIM) break;
	}

	tempSNo=ScanNum+1;
	vector<IsoCluster> XICList;
	size_t i, size=XICList_F.size();
	if(size>0)
	{
		for(i=size-1;i>0;i--) XICList.push_back(XICList_F[i]);
	}
	Interupt_times=0;
	while(tempSNo<Last)
	{
		tempSNo++;
		pt=getScanHeaderInfo(tempSNo);
		if(pt->m_data.msLevel!=1)  continue;
		pkl=getPeakList(tempSNo);
		if(pkl->getPeakCount()<=0) 
		{
			delete pt;
			continue;	
		}
		if(RefinePISO(it)) 
		{
			it.pTIC=GetMSTIC();
			it.RT=pt->getRetentionTimeSeconds()/60;
			XICList_F.push_back(it);
			Interupt_times=0;
		}
		else Interupt_times++;
		delete pkl;
		delete pt;
		if(Interupt_times>=IT_LIM) break;
	}

	vector<double> dMass;
	GetXICData(&XICList,&dMass,C);
	double m,s;
	CalMST(&dMass,&m,&s);
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s%d_%lf.%d.%d\n",sFName,dMass.size(),s,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",m*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",m);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);

	pkl=getPeakList(tempSNo);
	int PeakCount=pkl->getPeakCount();
	rampPeakInfoStruct *peakt;
	for(i=0;i<PeakCount;i++)
	{
		peakt=pkl->getPeak(i);
		if(peakt->intensity<1e-2) continue;
		fprintf(MGFfp,"%lf\t%lf\n",peakt->mz,peakt->intensity);
	}
	delete pkl;		
	fprintf(MGFfp,"END IONS\n");
	return true;
}

bool cRamp::SCalibrate(double C[5],FILE *MGFfp,long ScanNum)
{
	rampScanInfo *pt;
	if(ScanNum<First||ScanNum>Last) return 0;
	pt=getScanHeaderInfo(ScanNum);
	while(pt->m_data.msLevel!=2)  
	{
		delete pt;
		return 0;
	}
	delete pt;
	//refine the parent mz value in the MS scan
	bMass=pt->m_data.precursorMZ;
	charge=pt->m_data.precursorCharge;
	long tempSNo=ScanNum;	
	//seek to the last MS scan
	while(tempSNo>First)
	{
		tempSNo--;
		pt=getScanHeaderInfo(tempSNo);
		if(pt->m_data.msLevel==1)  break;
	}
	if(tempSNo<First) return false;
	pkl=getPeakList(tempSNo);
	//read the MS scan	
	if(pkl->getPeakCount()<=0) return false;
	//refine the pmz
	IsoCluster it;	
	if(!RefinePISO(it))
	{
		delete pkl;
		return false;
	}
	it.pTIC=GetMSTIC();
	it.RT=pt->getRetentionTimeSeconds()/60;
	delete pt;
	delete pkl;

	double mze=it.ISOP[0].dMass/MZE_SCALE;
	double tic=it.pTIC/TIC_SCALE;
	double rt=it.RT;
	double cmass=C[0]+C[1]*mze+C[2]*mze*mze+C[3]*mze*mze*tic+C[4]*rt;
	
	//output
	fprintf(MGFfp,"BEGIN IONS\nTITLE=%s.%d.%d\n",sFName,ScanNum,charge);
	//fprintf(MGFfp,"PMASS=%lf\n",cmass*charge-(charge-1)*ISO_DIFF[1]);//maybe wrong
	fprintf(MGFfp,"PMASS=%lf\n",cmass);//maybe wrong
	fprintf(MGFfp,"CHARGE=%d\n",charge);

	pkl=getPeakList(tempSNo);
	rampPeakInfoStruct *peakt;
	int PeakCount=pkl->getPeakCount();
	for(int i=0;i<PeakCount;i++)
	{
		peakt=pkl->getPeak(i);
		if(peakt->intensity<1e-2) continue;
		fprintf(MGFfp,"%lf\t%lf\n",peakt->mz,peakt->intensity);
	}
	delete pkl;		
	fprintf(MGFfp,"END IONS\n");
	return true;
}

int cRamp::Calibrate(double C[5],char *MGFfile)
{
	FILE *MGFfp;
	MGFfp=fopen(MGFfile,"w");
	if(MGFfp==NULL) return 0;
	int CT=0;
	for(long scannum=First;scannum<=Last;scannum++)
	{
		CT+=XICCalibrate(C,MGFfp,scannum);
	}
	fclose(MGFfp);
	return CT;
}

int cRamp::SimpleCalibrate(double C[5],char *MGFfile)
{
	FILE *MGFfp;
	MGFfp=fopen(MGFfile,"w");
	if(MGFfp==NULL) return 0;
	int CT=0;
	for(long scannum=First;scannum<=Last;scannum++)
	{
		CT+=SCalibrate(C,MGFfp,scannum);
	}
	fclose(MGFfp);
	return CT;
}

void cRamp::CalDefIsoDis(double mass,double IsoAreaT[6])
{
	IsoAreaT[0]=0.9937*exp(-0.0006693 *mass);

	double item=0.0006956*mass-0.03438;
	IsoAreaT[1]=item*exp(-item);

	item=0.0006899 *mass-0.018;
	IsoAreaT[2]=item*item*exp(-item)/2;

	item=0.0006673*mass+0.03829;
	IsoAreaT[3]=item*item*item*exp(-item)/6;
	
	item=0.0006462*mass+0.1133;
	IsoAreaT[4]=item*item*item*item*exp(-item)/24;

	item=0.0006522*mass+0.1314;
	IsoAreaT[5]=item*item*item*item*item*exp(-item)/120;
}


// END HENRY
