#include "StdAfx.h"
#include "mzXMLOut.h"
#include "OTrace.h"

using namespace std;

//void usage(const string& exename, const string& version) {
//	cout << endl << exename << " " << version << endl << endl;
//	cout << "Usage: " << exename << " [options] <raw file path> [<output file>]" << endl << endl
//		<< " Options\n"
//		<< "  --mzXML:         mzXML mode (default)" << endl
//		<< "  --mzML:          mzML mode (will use msconvert)" << endl
//		<< "      one of --mzXML or --mzML must be selected" << endl
//		<< endl
//		<< "  --centroid, -c: Centroid all scans (MS1 and MS2)" << endl
//		<< "      meaningful only if data was acquired in profile mode;"  << endl
//		<< "      default: off" << endl
//		<< "  [Advanced option, default OFF] --precursorFromFilterLine: only" << endl
//		<< "      try to get the precursor MZ value from the Thermo" << endl
//		<< "      \"filterline\" text; only use this if you have a good reason!" << endl
//		<< "      Otherwise, the program first will try to obtain a more accurate" << endl
//		<< "       mass from the \"Monoisotopic M/Z:\" \"trailer value\"" << endl 
//		<< "  --compress, -z: Use zlib for compressing peaks" << endl
//		<< "      default: off" << endl
//		<< "  --verbose, -v:   verbose" << endl
//		<< "  --gzip, -g:   gzip the output file (independent of peak compression)" << endl
//		<< endl
//		<< "  output file: (Optional) Filename for output file;" << endl
//		<< "      if not supplied, the output file will be created" << endl
//		<< "      in the same directory as the input file." << endl
//		<< endl
//		<< endl
//		<< "Example: convert input.raw file to output.mzXML, centroiding MS1 and MS2 scans" << endl << endl
//		<< "      " << exename << " --mzXML -c C:\\test\\input.raw c:\\test\\output.mzXML" << endl << endl
//		<< "Author: Natalie Tasman (SPC/ISB), with Jimmy Eng, Brian Pratt, and Matt Chambers," << endl
//		<< "      based on orignal work by Patrick Pedriolli." << endl;
//}


bool msconvert(ConverterArgs args,ThermoInterface thermoInterface)
{
	// force the program name to ReAdW,
	// regardless of what it was called on the command line
	const char *execName = "ReAdW"; 
	string version ="4.3.1";
	
	// try to init the thermo library
	/*if (!thermoInterface.initInterface()) 
	{
		printf("unable to interface with Thermo library\n");
		return -1;
	}

	if (!thermoInterface.setInputFile(args.inputFileName)) 
	{
		printf("unable to open %s with thermo interface", args.inputFileName);
		return -1;
	}*/

	thermoInterface.setCentroiding(args.centroidScans);
	thermoInterface.forcePrecursorFromFilter(args.forcePrecursorFromFilterLine);
	
	MassSpecXMLWriter* msWriter;
	if (args.mzMLMode)
	{
		msWriter = new mzMLWriter(execName, version, &thermoInterface);
	} 
	else if (args.mzXMLMode) 
	{
		msWriter = new mzXMLWriter(execName, version, &thermoInterface);
	}

	if (args.verbose) 
	{
		msWriter->setVerbose(true);
		thermoInterface.setVerbose(true);
	}

	if (!msWriter->setInputFile(args.inputFileName)) 
	{
		printf("unable to set input file %s\n ", args.inputFileName);
		return false;
	}

	if (!msWriter->setOutputFile(args.outputFileName)) 
	{
		printf("unable to open %s\n",args.outputFileName);
		return false;
	}

	msWriter->setCentroiding(args.centroidScans);
	msWriter->setCompression(args.compressScans);

	msWriter->writeDocument();
	
	if (args.verbose) 
	{
		int totalObsCharges = 0;
		printf("observed precursor charges");
		for (int i=0; i<int(thermoInterface.chargeCounts_.size()); i++) {
			printf("charge %d:%d\n", i,thermoInterface.chargeCounts_[i]);
			totalObsCharges += thermoInterface.chargeCounts_[i];
		}
		printf("total precursor charges: %d\n",totalObsCharges);
	}
	
	bool bRT=msWriter->IsOSuccess;
	delete msWriter;
	return bRT;
}

mzXMLOut::mzXMLOut(void)
{
}

mzXMLOut::~mzXMLOut(void)
{
}

