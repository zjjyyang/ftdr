// -*- mode: c++ -*-


/*
    File: MassSpecXMLWriter.cpp
    Description: instrument-independent common parent to mzXML and mzML writers.
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


#include "MassSpecXMLWriter.h"

#include <windows.h>
#include <atlconv.h>

#include <iostream>
#include <io.h>			// _findfirst, _open
#include <fcntl.h>		// _open 

#include "OTrace.h"

using namespace std;

std::string getComputerName()
{
	TCHAR computerName[MAX_COMPUTERNAME_LENGTH + 1];
	DWORD nameLen = (MAX_COMPUTERNAME_LENGTH+1) * sizeof(TCHAR);
	if (GetComputerName(computerName , &nameLen)) {
		USES_CONVERSION;
		return string(T2A(computerName));
	} else {
		return "";
	}
}

MassSpecXMLWriter::MassSpecXMLWriter(
									 const std::string& programName,
									 const std::string& programVersion,
									 InstrumentInterface* instrumentInterface)
	:   doCentroiding_(false),
		doDeisotoping_(false),
		doCompression_(false),
		shotgunFragmentation_(false),
		lockspray_(false),
		verbose_(false),
		doChecksum_(true)
{
	programName_ = programName;
	programVersion_ = programVersion;
	instrumentInterface_ = instrumentInterface;
	sha1Report_[0] = 0;
	// get hostname for the URIs
	hostName_ = getComputerName();
	printf("got computer name: %s\n",hostName_.c_str());
	//added by zhangjy on 2010.9.6
	IsOSuccess=false;
}

MassSpecXMLWriter::~MassSpecXMLWriter() {
	fout_.close();
}

void
MassSpecXMLWriter::setCentroiding(bool centroid) {
	doCentroiding_ = centroid;
	instrumentInterface_->setCentroiding(centroid);
}

void
MassSpecXMLWriter::setDeisotoping(bool deisotope) {
	doDeisotoping_ = deisotope;
	instrumentInterface_->setDeisotoping(deisotope);
}

void 
MassSpecXMLWriter::setCompression(bool compression) {
	doCompression_ = compression;
	instrumentInterface_->setCompression(compression);
}

void
MassSpecXMLWriter::setLockspray(bool ls) {
	lockspray_ = ls;
}

void
MassSpecXMLWriter::setShotgunFragmentation(bool sf) {
    shotgunFragmentation_ = sf;
}

void
MassSpecXMLWriter::setVerbose(bool verbose) {
	verbose_ = verbose;
}

bool
MassSpecXMLWriter::setInputFile(const std::string& inputFileName) {
	inputFileName_ = inputFileName;

	// TODO: deal with MassHunter data files correctly
	if (instrumentInterface_->instrumentInfo_.manufacturer_ == AGILENT) {
		// hack for trapper debugging
		inputFileNameList_.push_back(inputFileName+"AcqData\\Devices.xml");
		inputFileSHA1List_.push_back("dummy");
		return true;
	}

	//cout << "(got file URI: " << convertToURI(inputFileName_, hostName_) << ")" << endl;

	// run sha1

	if (instrumentInterface_->instrumentInfo_.manufacturer_ != WATERS) {
		// process one input file
		// sha1 input file
		if (verbose_ ) {
			printf("Calculating sha1-sum of %s\n",inputFileName.c_str());
		}
		sha1_.Reset();
		if ( !(sha1_.HashFile(inputFileName.c_str()))) {
			printf("Cannot open file %s  for sha-1 calculation\n",inputFileName);
			return false;// Cannot open file for sha1
		}
		sha1_.Final();
		sha1_.ReportHash(sha1Report_, SHA1::REPORT_HEX);
		if (verbose_) {
			printf("--done (sha1): %s\n", sha1Report_);
		}
	}
	else {
		// Waters
		// expect instrument interface has list of input filenames
		if (verbose_) {
			printf("Processing SHA1 for %d  input files.\n",inputFileNameList_.size());
		}
		for (vector<string>::size_type f=0; f<instrumentInterface_->inputFileNameList_.size(); ++f) {
			string inputFile = instrumentInterface_->inputFileNameList_[f];
			if (verbose_) {		
				printf("Calculating sha1-sum of %s\n",inputFile.c_str());
			}
			sha1Report_[0] = 0; // necessary-- fix SHA1 interface to be better
			sha1_.Reset();
			if ( !(sha1_.HashFile(inputFile.c_str()))) {
				printf("Cannot open file %S for sha-1 calculation\n",inputFile.c_str());
				return false;// Cannot open file for sha1
			}
			sha1_.Final();
			sha1_.ReportHash(sha1Report_, SHA1::REPORT_HEX);
			if (verbose_) {
				printf("--done (sha1): %s\n",sha1Report_);
			}

			// save these
			inputFileNameList_.push_back(inputFile);
			inputFileSHA1List_.push_back(sha1Report_);
		}

	}
	return true;
}

bool
MassSpecXMLWriter::setOutputFile(const std::string& outputFileName) {
	outputFileName_ = outputFileName;

	// Open output mzML file or die
	fout_.open(outputFileName.c_str());
	if( !fout_.good() ) {
		printf("Error opening output file %s\n",outputFileName.c_str());
		return false;
	}
	return true;
}


