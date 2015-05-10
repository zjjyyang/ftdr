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

#include "StdAfx.h"
#include "xsbrowsefolder.h"
#include <atlbase.h>

CXSBrowseFolder::CXSBrowseFolder(void)
:	m_style(0),
	m_root(NULL)
{

	ZeroMemory((PVOID)m_displayName, sizeof(m_displayName));
	ZeroMemory((PVOID)m_title, sizeof(m_title));

}

CXSBrowseFolder::~CXSBrowseFolder(void)
{
}

// Modifies the current style
DWORD CXSBrowseFolder::ModifyStyle(DWORD add, DWORD remove)
{
	m_style &= ~(remove);
	m_style |= add;

	return m_style;
}

// Returns the current style
DWORD CXSBrowseFolder::GetStyle(void)
{
	return m_style;
}

// Displays the dialog
CXSBrowseFolder::retCode CXSBrowseFolder::Show(HWND parent, LPTSTR pathBuffer)
{

	// Passed in a NULL pointer (!)
	ASSERT(pathBuffer);

	if (!pathBuffer)
		return RET_NOPATH;

	// Return value for the function
	retCode ret = RET_OK;

	LPITEMIDLIST itemIDList;

	// Set up the params
	BROWSEINFO browseInfo;

	// The dialog must have a parent
	ASSERT(parent);

	browseInfo.hwndOwner		= parent;
	browseInfo.pidlRoot			= m_root;
	browseInfo.pszDisplayName	= m_displayName;
	browseInfo.lpszTitle		= m_title;
	browseInfo.ulFlags			= m_style;
	browseInfo.lpfn				= NULL;
	browseInfo.lParam			= 0;
	
	// Show the dialog
	itemIDList = SHBrowseForFolder(&browseInfo);

	// Did user press cancel?
	if (!itemIDList)
		ret = RET_CANCEL;

	// Is everything so far?
	if (ret != RET_CANCEL) {
	
		// Get the path from the returned ITEMIDLIST
		if (!SHGetPathFromIDList(itemIDList, pathBuffer))
			ret = RET_NOPATH;
	
		// Now we need to free the ITEMIDLIST the shell allocated
		LPMALLOC	shellMalloc;
		HRESULT		hr;

		// Get pointer to the shell's malloc interface
		hr = SHGetMalloc(&shellMalloc);

		// Did it work?
		if (SUCCEEDED(hr)) {

			// Free the shell's memory
			shellMalloc->Free(itemIDList);

			// Release the interface
			shellMalloc->Release();

		}
	}
	
return ret;
}

// Set the title of the dialog
void CXSBrowseFolder::SetTitle(LPTSTR title)
{
	// NULL pointer (!)
	ASSERT(title);
	
	if (title)
		_tcscpy(m_title, title);
}
