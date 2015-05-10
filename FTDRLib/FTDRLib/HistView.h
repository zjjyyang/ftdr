#pragma once


// HistView �Ի���
#include "HistList.h"

class HistView : public CDialog
{
	DECLARE_DYNAMIC(HistView)

public:	
	HistView(CWnd* pParent = NULL);   // ��׼���캯��
	virtual ~HistView();
	bool DrawHist(CPaintDC *pDC,CRect &rect,myHistogram &Hist,bool Iscenter=false);
	void RemTail(char *buf);
	HistList HistsB;
	HistList HistsA;
	void scopy(char *buf,CString str);
	bool LoadData(char *fname);	
	bool SaveData(char *fname);	
	void SetData(HistList A[2]);	

// �Ի�������
	enum { IDD = IDD_DLGHIST };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnPaint();	
	afx_msg void OnBnClickedBtlast();
	afx_msg void OnBnClickedBtnext();
	afx_msg void OnBnClickedBtload();
	int m_Current;
	int m_total;
	CString m_rawfile;
	afx_msg void OnBnClickedBtsave();
};
