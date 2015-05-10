#pragma once

#include <string>
#include <process.h> // for system
#include "ThermoInterface.h"
#include "mzXMLWriter.h"
#include "mzMLWriter.h"
#include "ConverterArgs.h"
#include "MSUtilities.h"
#include "sysdepend.h"
#include "util.h"

class mzXMLOut
{
public:
	mzXMLOut(void);
	~mzXMLOut(void);	
};
bool msconvert(ConverterArgs args,ThermoInterface thermoInterface);