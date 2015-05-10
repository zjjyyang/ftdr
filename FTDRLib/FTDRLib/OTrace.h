// OTrace.h: interface for the OTrace class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_OTRACE_H__4EB75A68_410E_4D8D_8063_90B8960C21B1__INCLUDED_)
#define AFX_OTRACE_H__4EB75A68_410E_4D8D_8063_90B8960C21B1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifdef _DEBUG
#ifndef _OTRACE
#define _OTRACE
#endif //_OTRACE
#endif //_DEBUG

#ifdef	_OTRACE
#define OTRACE	_classOTrace.Print
#define OTRACE_TITLE	_classOTrace.SetTitle
#define GPOS _classOTrace.Getcurrent
#define SPOS _classOTrace.SetCurrent

	class OTrace  
	{ 
	public:
		OTrace();
		virtual ~OTrace();
		int Print(TCHAR *fmt, ...);
		void SetTitle(TCHAR *pTitle){SetConsoleTitle(pTitle);}			
		void Getcurrent(int *x,int *y);
		void SetCurrent(int row,int col);
	protected:
		HANDLE m_hStdOut;
	};
	extern OTrace _classOTrace;

#else //_OTRACE
	#define OTRACE
	#define OTRACE_TITLE
#endif //_OTRACE

#endif // !defined(AFX_OTRACE_H__4EB75A68_410E_4D8D_8063_90B8960C21B1__INCLUDED_)
