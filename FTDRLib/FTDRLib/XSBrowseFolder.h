/**************************************************************************************/
/* CXSBrowseFolder                                                                    */
/**************************************************************************************/
/* This is a simple class to wrap the SHBrowseForFolder function.                     */
/**************************************************************************************/
/* Written by Dana Holt, Xenos Software                                             */
/* http://www.xenossoftware.com                                                       */
/* This class is provided as-is, and carries no warranty or guarantee of any kind.    */
/* Use at your own risk.                                                              */
/**************************************************************************************/

#pragma once

class CXSBrowseFolder
{
public:

	enum retCode {

		RET_CANCEL = 0,
		RET_NOPATH,
		RET_OK
	};

public:
	CXSBrowseFolder(void);
	~CXSBrowseFolder(void);
protected:
	// Holds the current style
	DWORD m_style;
public:
	// Modifies the current style
	DWORD ModifyStyle(DWORD add, DWORD remove = 0);
	// Returns the current style
	DWORD GetStyle(void);
	// Displays the dialog
	CXSBrowseFolder::retCode Show(HWND parent, LPTSTR pathBuffer);
	// Set the title of the dialog
	void SetTitle(LPTSTR title);
protected:
	// Buffer to receieve the display name of the selected object
	TCHAR m_displayName[MAX_PATH];
	// Root item to start browsing at
	LPITEMIDLIST m_root;
	// Text to display above the tree view control
	TCHAR m_title[MAX_PATH];		
};
