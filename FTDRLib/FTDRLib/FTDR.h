// FTDR.h : main header file for the FTDR application
//

#if !defined(AFX_FTDR_H__290EFF06_6004_4EC4_926B_2C0E766FCF49__INCLUDED_)
#define AFX_FTDR_H__290EFF06_6004_4EC4_926B_2C0E766FCF49__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CFTDRApp:
// See FTDR.cpp for the implementation of this class
//

class CFTDRApp : public CWinApp
{
public:
	CFTDRApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFTDRApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CFTDRApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FTDR_H__290EFF06_6004_4EC4_926B_2C0E766FCF49__INCLUDED_)
