#if !defined(AFX_FILTERDLG_H__20074203_1EE8_4075_9575_586811A42A5D__INCLUDED_)
#define AFX_FILTERDLG_H__20074203_1EE8_4075_9575_586811A42A5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// FilterDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// FilterDlg dialog
#include "OutRecord.h"
#include "MResult.h"

class FilterDlg : public CDialog
{
// Construction
public:
	FilterDlg(CWnd* pParent = NULL);   // standard constructor
	SFilter f1;
	Filter f2;
	bool GetFilter(SFilter *sf,Filter *f);
	void SetFilter(SFilter *sf,Filter *f);
// Dialog Data
	//{{AFX_DATA(FilterDlg)
	enum { IDD = IDD_FILTER };
	BOOL	m_IsUseDefault;
	double	m_detCn;
	int		m_MaxMiss;
	int		m_MaxPepR;
	double	m_Met;
	double	m_MinPepScore;
	int		m_MaxRsp;
	double	m_XCorr1;
	double	m_XCorr2;
	double	m_XCorr3;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(FilterDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(FilterDlg)
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FILTERDLG_H__20074203_1EE8_4075_9575_586811A42A5D__INCLUDED_)
