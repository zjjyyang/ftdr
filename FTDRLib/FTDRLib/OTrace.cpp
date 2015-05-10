// OTrace.cpp: implementation of the OTrace class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "OTrace.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW

#ifndef _OTRACE
#define _OTRACE
#endif //_OTRACE

#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
#ifdef	_OTRACE
	OTrace::OTrace()
	{
		AllocConsole();
		SetConsoleTitle(_T("Work status"));
		m_hStdOut=GetStdHandle(STD_OUTPUT_HANDLE);
		COORD co={80,9999};
		SetConsoleScreenBufferSize(m_hStdOut,co);	
	}

	OTrace::~OTrace()
	{
		FreeConsole();
	}

	int OTrace::Print(TCHAR *fmt, ...)
	{
		TCHAR s[300];
		va_list argptr;
		int cnt;
		va_start(argptr,fmt);
		cnt=vsprintf(s,fmt,argptr);
		va_end(argptr);

		int slen=strlen(s);	
		DWORD cCharsWritten;
		if(m_hStdOut)
			WriteConsole(m_hStdOut,s,slen,&cCharsWritten,NULL);		
		return cnt;
	}
	void OTrace::Getcurrent(int *x,int *y)
	{
		CONSOLE_SCREEN_BUFFER_INFO ConsoleScreenBufferInfo;
		GetConsoleScreenBufferInfo(m_hStdOut,&ConsoleScreenBufferInfo);
		*x=ConsoleScreenBufferInfo.dwCursorPosition.Y;
		*y=ConsoleScreenBufferInfo.dwCursorPosition.X;
	}

	void OTrace::SetCurrent(int row,int col)
	{		
		COORD OutChar;//定义一个COORD类型的变量用于设置光标坐标
		OutChar.X = col;//光标设置为控制台屏幕开头
		OutChar.Y = row;	
		
		CONSOLE_SCREEN_BUFFER_INFO ConsoleScreenBufferInfo;
		GetConsoleScreenBufferInfo(m_hStdOut,&ConsoleScreenBufferInfo);
		int crow=ConsoleScreenBufferInfo.dwCursorPosition.Y;	
		int ccol=ConsoleScreenBufferInfo.dwCursorPosition.X;
		SetConsoleCursorPosition(m_hStdOut,OutChar);//设置光标位置
		//如果位置回退，刷新原来输出的字符
		DWORD cCharsWritten;
		int sizeB=(crow-row+1)*80;
		if(sizeB>0)	FillConsoleOutputCharacter(m_hStdOut,' ',sizeB,OutChar,&cCharsWritten);
		SetConsoleCursorPosition(m_hStdOut,OutChar);//设置光标位置		
	}

	OTrace _classOTrace;


#endif //_OTRACE
