# Introduction #
This is a practice build instruction for FTDR, all the steps have been tested.

FTDR build note:
Setting the libraries and include paths:
1. Configure the MSParser library.
Step 1: Download the MSParser from the Matrix Science Inc., the link is http://www.matrixscience.com/msparser.html, you must register your e-mail at first, the download link will be sent to your e-mail with a temporary user name and password.
The newest library version is for Visual studio 2008, both 64 bit version and 32 bit version, includes static the dynamic link libraries. The header files for include are also provided.
Step 2: Unzip the packages to a desired folder, such as “C:\Users\zjjyyang\DataAna\MSParser\MSParser32”
Step 3: setting the include path. In the IDE of visual studio 2008, open the FTDR project, open the project properties configuration, add the additional include path for MSParser in the C/C++ branch. In the General sub branch, the path can be select by browsing the Computer disk or input a string directly, in our setting, it should be “C:\Users\zjjyyang\DataAna\MSParser\MSParser32\vs2008\include”.
Step 4: setting the library path. In the project properties configuration panel, select the link branch, in the general sub branch, add the additional library dictionaries as “C:\Users\zjjyyang\DataAna\MSParser\MSParser32\vs2008\lib”.
Step 5: setting the linked library, for win 32 debug version, you should select MSParserD.lib, for the release version, MSParserS.lib or MSParserDS.lib may be better.
Note: you must select the sub folder in the MSParser64 for the win 64 program. Please try different configurations for the C runtime libraries to match the library selection in the C/C++ branch—Code generation—runtime library, and change the settings for MFC at the branch general for the whole project accordingly.
Most important, you cannot use the MSParser library in visual studio 2010, 2012 or 2013 as we know. I try it with many failures.
2. Configure the GSL library.
Step 1: download the GSL library for windows compile from
http://4fire.wordpress.com/2012/03/18/gsl-1-15-building-with-visual-studio-2010/
I recommend the 1.5 version for 2010, because I found some errors will present in the 1.6 version. Some problem about  the SSE2 instruction set will be appear if you do not have the Intel math libraries.
This website gives a step by step instruction, and anyone can do it easily. In our setting, the build  library is in the path “C:\Users\zjjyyang\DataAna\GSL\gsl-1.15-vc10”.
Step 2. Configuration for FTDR. Add include path “C:\Users\zjjyyang\DataAna\GSL\gsl-1.15-vc10”, add library path “C:\Users\zjjyyang\DataAna\GSL\gsl-1.15-vc10\build.vc10\lib\Win32\Release”. Note, we use the static libraries here.
Step 3. Add the linked library: gsl.lib and cblas.lib.
3. Configure the boost library.
Step 1： download the correct version of boost library. Now, the library has been complied for different visual studio version, both 32 and 64 bit version are available. The website is http://www.boost.org/users/history/version_1_55_0.html.
Please select the necessary version from “Other Downloads”, you will be redirect to the SourceForge.
Here, we download the 32 bit version for VC 9.
Step 2: install the boost library, just click the setup file, the default folder is c:\local\boost\_1\_55\_0.
Step 3: configuration in the VS 2008, include path: C:\local\boost\_1\_55\_0, library path: C:\local\boost\_1\_55\_0\lib32-msvc-9.0.
Step 4: add the correct library you used in the program, here, we do not specify, let it select by the system.
4. Configure the ProteoWizard Library.
Step 1: download ProteoWizard from http://proteowizard.sourceforge.net/
Several options are provided; we select the source MSVC build (no 3th party reader support) because FTDR only used the mzXML and mzML library.
Step 2: unzip this library, copy necessary **.hpp and**.cpp files for the mzML and mzXML IO, we do not compile the library in the FTDR project, but use the source code directly. All the configurations have been tested in this project, it works well. Of course, the user can also try to use the static or dll libraries.
5. Download the zlib library: http://gnuwin32.sourceforge.net/packages/zlib.htm
Just install the binaries.
6. Compile FTDR.