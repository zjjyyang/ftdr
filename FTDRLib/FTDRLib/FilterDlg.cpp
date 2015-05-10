// FilterDlg.cpp : implementation file
//

#include "stdafx.h"
#include "FTDR.h"
#include "FilterDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// FilterDlg dialog


FilterDlg::FilterDlg(CWnd* pParent /*=NULL*/)
	: CDialog(FilterDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(FilterDlg)
	m_IsUseDefault = TRUE;
	m_detCn = 0.1f;
	m_MaxMiss = 2;
	m_MaxPepR = 10;
	m_Met = 15.0f;
	m_MinPepScore = 0.0f;
	m_MaxRsp = 4;
	m_XCorr1 = 1.5f;
	m_XCorr2 = 2.5f;
	m_XCorr3 = 3.0f;
	//}}AFX_DATA_INIT
}


void FilterDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(FilterDlg)
	DDX_Check(pDX, IDC_CISUDPSF, m_IsUseDefault);
	DDX_Text(pDX, IDC_EDDETCN, m_detCn);
	DDX_Text(pDX, IDC_EDMAXMISS, m_MaxMiss);
	DDX_Text(pDX, IDC_EDMAXPEPR, m_MaxPepR);
	DDX_Text(pDX, IDC_EDMET, m_Met);
	DDX_Text(pDX, IDC_EDMINPEPS, m_MinPepScore);
	DDX_Text(pDX, IDC_EDRSP, m_MaxRsp);
	DDX_Text(pDX, IDC_EDXCORR1, m_XCorr1);
	DDX_Text(pDX, IDC_EDXCORR2, m_XCorr2);
	DDX_Text(pDX, IDC_EDXCORR3, m_XCorr3);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(FilterDlg, CDialog)
	//{{AFX_MSG_MAP(FilterDlg)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// FilterDlg message handlers

void FilterDlg::OnOK() 
{
	// TODO: Add extra validation here
	UpdateData();
	f1.XCorr[0]=m_XCorr1;
	f1.XCorr[1]=m_XCorr2;
	f1.XCorr[2]=m_XCorr3;
	f1.detCn=m_detCn;
	f1.RSp=m_MaxRsp;
	f2.MinPepScore=m_MinPepScore;
	f2.IsUseDefault=m_IsUseDefault;
	f2.MET=m_Met*1e-6f;
	f2.MaxRank=m_MaxPepR;
	f2.MaxMiss=m_MaxMiss;	
	CDialog::OnOK();
}

bool FilterDlg::GetFilter(SFilter *sf,Filter *f)
{
	*sf=f1;
	*f=f2;
	return true;
}

void FilterDlg::SetFilter(SFilter *sf,Filter *f)
{
	f1=*sf;
	f2=*f;
}

BOOL FilterDlg::OnInitDialog() 
{
	CDialog::OnInitDialog();
	
	// TODO: Add extra initialization here
	m_XCorr1=f1.XCorr[0];
	m_XCorr2=f1.XCorr[1];
	m_XCorr3=f1.XCorr[2];
	m_detCn=f1.detCn;
	m_MaxRsp=f1.RSp;
	m_MinPepScore=f2.MinPepScore;
	m_IsUseDefault=f2.IsUseDefault;
	m_Met=f2.MET*1e6f;
	m_MaxPepR=f2.MaxRank;
	m_MaxMiss=f2.MaxMiss;	
	UpdateData(false);
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}
