// mzXML.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "mzXML.h"

#include "cramp.hpp"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 唯一的应用程序对象

CWinApp theApp;

using namespace std;

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;

	// 初始化 MFC 并在失败时显示错误
	if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	{
		// TODO: 更改错误代码以符合您的需要
		_tprintf(_T("错误: MFC 初始化失败\n"));
		nRetCode = 1;
	}
	else
	{
		// TODO: 在此处为应用程序的行为编写代码。
		cRamp mzXMLFile(argv[1]);
		int Last=mzXMLFile.getLastScan();
		int First=1;
		rampScanInfo *pt;
		for(int i=First;i<=Last;i++)
		{
			pt=mzXMLFile.getScanHeaderInfo(i);
			if(pt->m_data.msLevel==1) cout<<pt->m_data.peaksCount<<"\t"<<pt->m_data.scanType<<endl;
		}
	}

	return nRetCode;
}
