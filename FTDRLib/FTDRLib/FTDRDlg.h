// FTDRDlg.h : header file
//

#if !defined(AFX_FTDRDLG_H__7858DA13_D3C4_4729_9076_230D9B370B3D__INCLUDED_)
#define AFX_FTDRDLG_H__7858DA13_D3C4_4729_9076_230D9B370B3D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CFTDRDlg dialog
#include "OTrace.h"
#include "xsbrowsefolder.h"
#include "alQuant.h"
//#include "afxwin.h"
#include <vector>
#include "afxwin.h"
#include "afxcmn.h"
using std::vector;

class CFTDRDlg : public CDialog
{
// Construction
public:
	CString ErrorMessage[8];
	bool Write_Error_code(CReport *rt,CStdioFile *dFile);
	void Add_Error_Msg(CReport *rt);
	void AddRawItem(CString raw);
	void AddMsg(CString msg);
	void WriteReport();	
	void WriteReport(int i);
	void WriteReport(CReport &rt);
	void ADRList();	
	void AddRitem(CString str,CString Ext);
	void InitialCtrl();
	//void AdjustSVMFName();
	CFTDRDlg(CWnd* pParent = NULL);	// standard constructor
	LPTSTR GetNTS(CString cString);
	//OTrace dbg;
	CString m_ResultFs,m_RawFs;
	ThreadParm *Parm;
	void PackParm();
	int ValidNum,ThreadNum;
	int CurrentFP;
	void GetFName(CString Path, CString &fn);
	CString GetFileExt(CString fname);
	CString GetPureFName(CString fname);
	int detProg,prog;
	vector<CReport> RepList;
	//int MsgLine;
	CString m_sRawFile;
	//CString svm_file;

	//to process different type messages
	void MsgThreadEnd(LPARAM lParam);
	bool MsgTaskFinshed(LPARAM lParam);
	void MsgRawNotices(WPARAM wParam,LPARAM lParam);

	//advance paremeters setting
	AdvPars APars;
	void UpdateEffective();

	//save experiment designs
	bool SaveEXPtoFile(CString &SFName);
	bool LoadEXPFromFile(CString &SFName);

// Dialog Data
	//{{AFX_DATA(CFTDRDlg)
	enum { IDD = IDD_FTDR_DIALOG };
	CListCtrl	m_RawList;
	//CListBox	m_MsgList;
	CComboBox	m_ThreadNum;
	CProgressCtrl	m_WorkPro;
	CListCtrl	m_RList;
	CString	m_OutPath;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFTDRDlg)
	public:
	virtual BOOL DestroyWindow();
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CFTDRDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnCal();
	afx_msg void OnBtAddFile();
	afx_msg void OnBtAddFiles();
	afx_msg void OnBtAddPath();
	afx_msg void OnBtAddPathes();
		afx_msg void OnBtRemOne();
	afx_msg void OnBtRemAlll();
	afx_msg void OnBtoutpath();
	afx_msg void OnBtfilter();
	afx_msg LRESULT OnDisplay(WPARAM wParam,LPARAM lParam);
	afx_msg void OnBtaddsraw();
	afx_msg void OnBtbatchraw();
	afx_msg void OnBtremsraw();
	afx_msg void OnBtremaraw();
	afx_msg void OnBnClickedBtviewhist();
	afx_msg void OnBtRemAll() ;
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:	
	CButton m_RPO;
	CButton m_OPD;	
public:
	//afx_msg void OnBnClickedBtModels();
public:
	CComboBox MS1OFormat;
	afx_msg void OnBnClickedBtaddsetting();
	afx_msg void OnBnClickedBtloadds();
	afx_msg void OnBnClickedBtsaveds();
	CComboBox ModelSel;
	afx_msg void OnCbnSelchangeCombmodel();
	CButton StatusPar;
	CComboBox MS2Type;
	CRichEditCtrl MSGPostBox;
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FTDRDLG_H__7858DA13_D3C4_4729_9076_230D9B370B3D__INCLUDED_)
