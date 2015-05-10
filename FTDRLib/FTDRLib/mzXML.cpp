// mzXML.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "mzXML.h"

#include "cramp.hpp"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// Ψһ��Ӧ�ó������

CWinApp theApp;

using namespace std;

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;

	// ��ʼ�� MFC ����ʧ��ʱ��ʾ����
	if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	{
		// TODO: ���Ĵ�������Է���������Ҫ
		_tprintf(_T("����: MFC ��ʼ��ʧ��\n"));
		nRetCode = 1;
	}
	else
	{
		// TODO: �ڴ˴�ΪӦ�ó������Ϊ��д���롣
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
