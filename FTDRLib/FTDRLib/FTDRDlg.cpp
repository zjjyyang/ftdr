// FTDRDlg.cpp : implementation file
//

#include "stdafx.h"
#include "FTDR.h"
#include "FTDRDlg.h"
#include "FilterDlg.h"
#include "WorkThread.h"
#include <gsl/gsl_multifit.h>
#include "HistView.h"
#include "AdvSetDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFTDRDlg dialog

CFTDRDlg::CFTDRDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CFTDRDlg::IDD, pParent)	
{
	//{{AFX_DATA_INIT(CFTDRDlg)
	m_OutPath = _T("Please Select the output Path!");
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CFTDRDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CFTDRDlg)
	DDX_Control(pDX, IDC_LIST4, m_RawList);
	//DDX_Control(pDX, IDC_LIST3, m_MsgList);
	DDX_Control(pDX, IDC_CBTHREAD, m_ThreadNum);
	DDX_Control(pDX, IDC_PROGRESS1, m_WorkPro);
	DDX_Control(pDX, IDC_LIST1, m_RList);
	DDX_Text(pDX, IDC_ETOUTPATH, m_OutPath);	
	//}}AFX_DATA_MAP
	DDX_Control(pDX, IDC_CHECK3, m_RPO);
	DDX_Control(pDX, IDC_CHECK4, m_OPD);
	DDX_Control(pDX, IDC_COMBO1, MS1OFormat);
	DDX_Control(pDX, IDC_COMBMODEL, ModelSel);
	DDX_Control(pDX, IDC_CKSTSPAR, StatusPar);
	DDX_Control(pDX, IDC_COMMS2, MS2Type);
	DDX_Control(pDX, IDC_RICHEDITMSG, MSGPostBox);
}

BEGIN_MESSAGE_MAP(CFTDRDlg, CDialog)
	//{{AFX_MSG_MAP(CFTDRDlg)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_CAL, OnCal)
	ON_BN_CLICKED(IDC_BTADDONEMASCOT, OnBtAddFile)
	ON_BN_CLICKED(IDC_BTBATCHMASCOT, OnBtAddFiles)
	ON_BN_CLICKED(IDC_BTSSEQUEST, OnBtAddPath)
	ON_BN_CLICKED(IDC_BTBSEQUEST, OnBtAddPathes)	
	ON_BN_CLICKED(IDC_BTREMONE, OnBtRemOne)
	ON_BN_CLICKED(IDC_BTREMALL, OnBtRemAll)	

	ON_BN_CLICKED(IDC_BTOUTPATH, OnBtoutpath)
	ON_BN_CLICKED(IDC_BTFILTER, OnBtfilter)
	ON_MESSAGE(WM_DISPLAY,OnDisplay)
	ON_BN_CLICKED(IDC_BTADDSRAW, OnBtaddsraw)
	ON_BN_CLICKED(IDC_BTBATCHRAW, OnBtbatchraw)
	ON_BN_CLICKED(IDC_BTREMSRAW, OnBtremsraw)
	ON_BN_CLICKED(IDC_BTREMARAW, OnBtremaraw)
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BTVIEWHIST, &CFTDRDlg::OnBnClickedBtviewhist)
	//ON_BN_CLICKED(IDC_BT_MODELS, &CFTDRDlg::OnBnClickedBtModels)
	ON_BN_CLICKED(IDC_BTADDSETTING, &CFTDRDlg::OnBnClickedBtaddsetting)
	ON_BN_CLICKED(IDC_BTLOADDS, &CFTDRDlg::OnBnClickedBtloadds)
	ON_BN_CLICKED(IDC_BTSAVEDS, &CFTDRDlg::OnBnClickedBtsaveds)
	ON_CBN_SELCHANGE(IDC_COMBMODEL, &CFTDRDlg::OnCbnSelchangeCombmodel)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFTDRDlg message handlers

BOOL CFTDRDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.
	SetWindowPos(&CWnd::wndNoTopMost, 0, 10, 0, 70, SWP_NOSIZE);

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	
	// TODO: Add extra initialization here
	InitialCtrl();
	Parm=new ThreadParm;
	detProg=0;
	prog=0;
//	MsgLine=0;
	//svm_file="svm.txt";
	//////////////////////////////////////////
	ErrorMessage[0]="open raw file failure!";
	ErrorMessage[1]="no records for model building loaded,possilbe reasons:result file format (such as mascot html files) is not supported!";
	ErrorMessage[2]="create output file failure!";
	ErrorMessage[3]="calibrating model building failure, too less validated identifications for model building!";
	ErrorMessage[4]="can notopen the results file or path no out files!";
	ErrorMessage[5]="Incorrect database search result format, only support sequest path, mascot dat/html, pepxml format.";
	ErrorMessage[6]="Not output global ms1 as mzML or mzXML, failure or not select.";
	ErrorMessage[7]="Not support model, please check the model sleection combox.";
	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CFTDRDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CFTDRDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CFTDRDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CFTDRDlg::PackParm()
{
	UpdateData();
	ThreadNum=m_ThreadNum.GetCurSel()+1;	
	Parm->IsEnd=false;
	Parm->IsOpCData=m_OPD.GetCheck();
	Parm->IsOPSVMModel=m_RPO.GetCheck();
	Parm->IsUseStatsPar=StatusPar.GetCheck();
	//if(Parm->IsOPSVMModel) Parm->SVMModelFile=svm_file;
	int sel=MS1OFormat.GetCurSel();
	if(sel==0)Parm->ms1format=MS1_NON;
	else if(sel==1)Parm->ms1format=MS1_mzXML;	
	else if(sel==2)Parm->ms1format=MS1_mzML;
	
	sel=MS2Type.GetCurSel();
	if(sel==0)Parm->ms2format=MS2NON;
	else if(sel==1)Parm->ms2format=MS2MGF;	
	else if(sel==2)Parm->ms2format=MS2EMGF;	

	sel=ModelSel.GetCurSel();
	if(sel==0) Parm->ModelType=MODEL_NON;
	else if(sel==1) Parm->ModelType=MODEL_LIN;
	else if(sel==2) Parm->ModelType=MODEL_SVM;
	else Parm->ModelType=MODEL_ERROR;

	Parm->hEventP=CreateEvent(NULL, TRUE, FALSE, NULL);
	Parm->hEventE=CreateEvent(NULL, TRUE, TRUE, NULL);
	Parm->hEventB=CreateEvent(NULL, TRUE, TRUE, NULL);
	int RN=m_RList.GetItemCount();	
	ValidNum=RN;
	int i=m_RawList.GetItemCount();
	ValidNum=ValidNum>i?i:ValidNum;
	if(ValidNum<RN)
	{
		AddMsg("Number of result files is more than raw files, the tail will be trunked!");
	}
	
	if(ValidNum<i)
	{
		AddMsg("Number of Raw files is more than result files, the tail will be trunked!");
	}
	if(ThreadNum>ValidNum) ThreadNum=ValidNum;
	Parm->AllocRepBuff(ThreadNum);
	CurrentFP=0;
}
//
//void CFTDRDlg::AdjustSVMFName()
//{
//	CString item;
//	item.Format("_%d",CurrentFP);
//
//	int si=svm_file.GetLength()-1;
//	int idx=svm_file.ReverseFind('.');
//	if(idx==-1) Parm->SVMModelFile=svm_file+item;
//	else
//	{
//		CString tmp=svm_file.Mid(idx);
//		Parm->SVMModelFile=svm_file.Mid(0,idx);
//		Parm->SVMModelFile+=item;
//		Parm->SVMModelFile+=tmp;
//	}
//}

void CFTDRDlg::GetFName(CString Path, CString &fn)
{
	int idx=Path.GetLength()-1;
	while(idx>=0&&Path[idx]!='\\') idx--;
	idx++;
	fn=Path.Mid(idx);
}

void CFTDRDlg::OnCal() 
{
	// TODO: Add your control notification handler code here
	UpdateData();
	RepList.clear();
	if(m_OutPath=="Please Select the output Path!")
	{
		MessageBox("Please input your output path!");
		return;
	}	
	if(m_RawList.GetItemCount()<=0)
	{
		MessageBox("No Raw file is provided!");
		return;
	}

	if(m_RList.GetItemCount()<=0)
	{
		MessageBox("No database search result file is provided!");
		return;
	}
	PackParm();
	int i;
	for(i=0;i<ThreadNum;i++)
	{
		WaitForSingleObject(Parm->hEventB,INFINITE);
		Parm->ThreadID=i;
		ResetEvent(Parm->hEventB);
		AfxBeginThread(OutSumThread,Parm);		
	}
	CString st1,st3,st4;
	int issuccess=SetEvent(Parm->hEventE);
	CString rtype,dtype;
	int type;
	CString sMsg;

	for(i=0;i<ValidNum;i++)//process each thread
	{
		if(CurrentFP>=ThreadNum) break;
		//AdjustSVMFName();
		st1=m_RList.GetItemText(i,1);	
		st4=m_RawList.GetItemText(i,1);
		rtype=m_RList.GetItemText(i,2);		
		st3=m_OutPath;
		if((rtype=="htm"||rtype=="dat"))type=R_TYPE_MASCOT;//Mascot html or dat	
		else if(rtype=="path")type=R_TYPE_SEQUEST;//Sequest output path
		else if(rtype=="pepxml") type=R_TYPE_PEPXML;//PepXML
		else if(rtype=="massmatrix") type=R_TYPE_MASSMATRIX;
		else
		{
			OTRACE("result type for %d row is not support!\n", i);
			sMsg.Format("result type for %d row is not support!", i);
			AddMsg(sMsg);
			continue;
		}
		WaitForSingleObject(Parm->hEventE,INFINITE);
		OTRACE("Sending %d: %s!\n", CurrentFP,st1);		
		sMsg.Format("Sending %d: %s!\n", CurrentFP,(LPCTSTR)st1);
		AddMsg(sMsg);

        CurrentFP++;		
		Parm->GL_CriticalSection.Lock();	
		Parm->ThreadID=CurrentFP-1;
		Parm->RFile=st1;	//assign the result file or path
		Parm->RawFile=st4;	//assign the raw file
		Parm->sOutPath=st3;	//assign the outpath
		Parm->type=type;   //assign the result type
		Parm->IsDeal=false;//initial the process status
		SetEvent(Parm->hEventP); //Synchronization control
		ResetEvent(Parm->hEventE);
		Parm->GL_CriticalSection.Unlock();
	}
	((CButton *)GetDlgItem(IDC_CAL))->EnableWindow(false);
	detProg=100/ValidNum; //initial the progress bar
	prog=0;	
	m_WorkPro.SetPos(prog);
}

LPTSTR CFTDRDlg::GetNTS(CString cString)
{
   LPTSTR lpsz = new TCHAR[cString.GetLength()+1];
   _tcscpy(lpsz, cString);
   return lpsz;
}

CString CFTDRDlg::GetFileExt(CString fname)
{
	CString Ext;
	int idx=fname.ReverseFind('.');
	if(idx!=-1)
	{
		Ext=fname.Mid(idx+1);
		if(Ext=="xml")	Ext="pepxml";
		else if(Ext=="csv") Ext="massmatrix";
	}
	return Ext;
}

CString CFTDRDlg::GetPureFName(CString fname)
{
	CString Ext;
	int idx=fname.ReverseFind('.');
	if(idx!=-1)
	{
		Ext=fname.Mid(0,idx);
	}
	idx=Ext.ReverseFind('\\');
	if(idx!=-1)
	{
		Ext=Ext.Mid(idx+1);
	}
	return Ext;
}

void CFTDRDlg::OnBtAddFile() 
{
	// TODO: Add your control notification handler code here
	TCHAR szPathBuffer[_MAX_PATH];
	TCHAR szDrive[_MAX_DRIVE];
	TCHAR szDir[_MAX_DIR];

	CString sFileFilter(_T("Database search result Files (Mascot or PepXML) (*.htm;*.dat;*.pepXML;*.xml;*.csv) | *.htm;*.dat;*.pepXML;*.xml;*.csv||"));
	_tsplitpath(m_ResultFs, szDrive, szDir, NULL, NULL);
	_tmakepath(szPathBuffer, szDrive, szDir, NULL, NULL);

	CFileDialog dlg(TRUE, 
					_T("htm"), 
					_T("*.htm;*.dat;*.pepXML;*.xml;*.csv"), 
					OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_EXTENSIONDIFFERENT | OFN_EXPLORER,
					sFileFilter,
					this );
	dlg.m_ofn.lpstrTitle = _T("Open database search result File");
	dlg.m_ofn.lpstrDefExt = _T("htm");
	dlg.m_ofn.lpstrInitialDir = szPathBuffer;
    
	if(dlg.DoModal()==IDOK)
	{
         m_ResultFs=dlg.GetPathName();
		 CString sExt;		 
		 if(m_ResultFs.Find("pep.xml")==m_ResultFs.GetLength()-7) sExt="pepxml";
		 else 
		 {
			 sExt=dlg.GetFileExt();			
			 sExt.MakeLower(); 
			 if(sExt=="csv") sExt="massmatrix";
		 }
		 //if(sExt=="dat"||sExt=="htm"||sExt=="html"||sExt=="pepxml")) 
		 AddRitem(m_ResultFs,sExt);		
		 //else MessageBox("Incorrect file format or extentions!");
	}
}

void CFTDRDlg::OnBtAddFiles() 
{
	// TODO: Add your control notification handler code here
   CXSBrowseFolder browse;
   browse.SetTitle("Open the directory with *.htm, *.html, *.dat, *.pepxml, *.csv or *.pepXML files");
   char temppath[MAX_PATH];
   CString ExtENUM[6]={
	   "*.htm",
	   "*.html",
	   //"*.HTM",
	   //"*.HTML",
	  // "*.pepXML",
	   "*.pepxml",
	   "*.pep.xml",
	  // "*.PEPXML",
	   "*.dat",
	   //"*.DAT"
	   "*.csv",
   };
   browse.ModifyStyle(BIF_EDITBOX|BIF_RETURNONLYFSDIRS);  
   if(browse.Show(this->GetSafeHwnd(),temppath)==2)
   {
	   for(int i=0;i<5;i++)
	   {
		  HANDLE hFind;
		  WIN32_FIND_DATA fd;
		  CString strFileSpec =temppath;
		  if (strFileSpec.Right (1) != "\\") strFileSpec += "\\";
		  strFileSpec +=ExtENUM[i] ;
		  if ((hFind =::FindFirstFile ((LPCTSTR) strFileSpec, &fd))==INVALID_HANDLE_VALUE)  continue;
		  CWaitCursor wait;
		  do 
		  {
			  CString strFileName = (LPCTSTR) &fd.cFileName;
			  CString temp=temppath;
			  temp+="\\";
			  temp=temp+GetNTS(fd.cFileName);
			  temp.Replace(":\\\\",":\\");
			  CString Ext= GetFileExt(temp);		
			  AddRitem(temp,Ext);		
		  } while (::FindNextFile (hFind, &fd));
		  ::FindClose (hFind);	
	   }	 
   }   
}

void CFTDRDlg::OnBtAddPath() 
{
	// TODO: Add your control notification handler code here
   CXSBrowseFolder browse;
   browse.SetTitle("Open the directory with *.out files");
   char temppath[MAX_PATH];
   browse.ModifyStyle(BIF_EDITBOX|BIF_RETURNONLYFSDIRS);  
   if(browse.Show(this->GetSafeHwnd(),temppath)==2)
   {
	   CString stemp=temppath;
	   AddRitem(stemp,"path");
   }
}

void CFTDRDlg::OnBtAddPathes() 
{
	// TODO: Add your control notification handler code here
   CXSBrowseFolder browse;
   browse.SetTitle("select sub-directories in this directory.");
   char temppath[MAX_PATH];
   browse.ModifyStyle(BIF_EDITBOX|BIF_RETURNONLYFSDIRS);  
   if(browse.Show(this->GetSafeHwnd(),temppath)==2)
   {
	  CString temp;
	  HANDLE hFind;
	  WIN32_FIND_DATA fd;
	  CString strFileSpec =temppath;
	  if (strFileSpec.Right (1) != "\\") strFileSpec += "\\";
	  strFileSpec += "*.*";
	  if ((hFind =::FindFirstFile ((LPCTSTR) strFileSpec, &fd))==INVALID_HANDLE_VALUE) 
	  { 
		  return ;
	  }
	  CWaitCursor wait;
	  do 
	  {
		  if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) 
		  {
			  CString strFileName = (LPCTSTR) &fd.cFileName;
			  if ((strFileName != ".") && (strFileName != "..")&& (fd.dwFileAttributes != 22))
			  {
				  temp=temppath;
				  temp+="\\";
				  temp=temp+GetNTS(fd.cFileName);
				  temp.Replace(":\\\\",":\\");
				  AddRitem(temp,"path");				  			
			  }
		  }
	  } while (::FindNextFile (hFind, &fd));
	  ::FindClose (hFind);	  
	  ////////////////////sort path 	 
   } 
}

void CFTDRDlg::OnBtRemOne() 
{
	// TODO: Add your control notification handler code here
	CString *pList;
	CString *tList;
	int count=m_RList.GetItemCount();
	if(count<=0) return;
	pList=new CString[count];
	tList=new CString[count];
	int i;
	for(i=0;i<count;i++)
	{		
		pList[i]=m_RList.GetItemText(i,1);
		tList[i]=m_RList.GetItemText(i,2);
	}
	POSITION pos = m_RList.GetFirstSelectedItemPosition();
	if (pos == NULL) return;		
	else
	{
		while (pos)
		{
			int nItem = m_RList.GetNextSelectedItem(pos);
			pList[nItem].Empty();
			//m_RList.DeleteItem(nItem);
			//OTRACE("Item %d was selected!\n", nItem);		
		}
	}
	m_RList.DeleteAllItems();
	for(i=0;i<count;i++)
	{
		if(pList[i].IsEmpty()) continue;	
		AddRitem(pList[i],tList[i]);
	}
	delete []pList;
	delete []tList;
	//ADRList();
}

void CFTDRDlg::OnBtRemAll() 
{
	// TODO: Add your control notification handler code here
	m_RList.DeleteAllItems();	
}

void CFTDRDlg::OnBtoutpath() 
{
	// TODO: Add your control notification handler code here
   UpdateData();
   CXSBrowseFolder browse;
   browse.SetTitle("Open the directory with *.htm files");
   char temppath[MAX_PATH];
   browse.ModifyStyle(BIF_EDITBOX|BIF_RETURNONLYFSDIRS);  
   if(browse.Show(this->GetSafeHwnd(),temppath)==2)
   {
	   m_OutPath=temppath;
	   UpdateData(false);
   }	
}

void CFTDRDlg::OnBtfilter() 
{
	// TODO: Add your control notification handler code here
	FilterDlg dlg;
	dlg.SetFilter(&Parm->dFilter1,&Parm->dFilter);
	if(dlg.DoModal()==IDOK)
	{
		dlg.GetFilter(&Parm->dFilter1,&Parm->dFilter);
	}	
}

void CFTDRDlg::InitialCtrl()
{
	LV_COLUMN lvc;
	lvc.mask	= LVCF_FMT | LVCF_TEXT | LVCF_WIDTH;
	lvc.fmt		= LVCFMT_LEFT;
	lvc.cx		= 40;

	CString sHeading = _T("No.");
	lvc.pszText		= (LPTSTR)(LPCTSTR)sHeading;
	lvc.cchTextMax	= sHeading.GetLength();

	m_RList.InsertColumn(0,&lvc);
	m_RawList.InsertColumn(0,&lvc);

	sHeading = _T("File or Path");
	lvc.cx		= 310;
	lvc.pszText		= (LPTSTR)(LPCTSTR)sHeading;
	lvc.cchTextMax	= sHeading.GetLength();

	m_RList.InsertColumn(1,&lvc);	
	lvc.cx		= 370;
	sHeading = _T("File");
	lvc.pszText		= (LPTSTR)(LPCTSTR)sHeading;
	lvc.cchTextMax	= sHeading.GetLength();
	m_RawList.InsertColumn(1,&lvc);

	sHeading = _T("Type");
	lvc.cx		= 60;
	lvc.pszText		= (LPTSTR)(LPCTSTR)sHeading;
	lvc.cchTextMax	= sHeading.GetLength();

	m_RList.InsertColumn(2,&lvc);
	/////////////////////////////////
	for(int i=1;i<11;i++)
	{
		sHeading.Format("%d",i);
		m_ThreadNum.AddString(sHeading);
	}
	m_ThreadNum.SetCurSel(0);
	m_WorkPro.SetRange(0,100);
	m_WorkPro.SetPos(100);	

	MS1OFormat.AddString("Not output");
	MS1OFormat.AddString("mzXML");
	MS1OFormat.AddString("mzML");
	MS1OFormat.SetCurSel(0);

	MS2Type.AddString("Non");
	MS2Type.AddString("mgf");
	MS2Type.AddString("Extend mgf");
	MS2Type.SetCurSel(1);

	m_RList.SetExtendedStyle(m_RList.GetExtendedStyle()|LVS_EX_FULLROWSELECT); 
	m_RawList.SetExtendedStyle(m_RawList.GetExtendedStyle()|LVS_EX_FULLROWSELECT);

	ModelSel.AddString("Non");
	ModelSel.AddString("Linear");
	ModelSel.AddString("svm");
	ModelSel.SetCurSel(2);

	StatusPar.SetCheck(true);
	m_OPD.SetCheck(true);
	m_RPO.SetCheck(true);
}

void CFTDRDlg::AddRitem(CString str,CString Ext)
{
	int i=m_RList.GetItemCount();	
	LV_ITEM lvi;
	lvi.mask= LVIF_TEXT; 
	CString sTemp;
	lvi.iItem		= i;
	lvi.iSubItem	= 0;
	sTemp.Format("%d",i+1);
	lvi.pszText		= (LPTSTR)(LPCTSTR)sTemp;
	lvi.cchTextMax	= sTemp.GetLength();
	m_RList.InsertItem(&lvi); 

	lvi.iSubItem	= 1;
	sTemp=str;
	lvi.pszText		= (LPTSTR)(LPCTSTR)sTemp;
	lvi.cchTextMax	= sTemp.GetLength(); 
	m_RList.SetItem(&lvi);
	

	//if(type==0) sTemp="htm";
	//else if(type==1) sTemp="Path";
	//else if(type==2) sTemp="dat";
	//else if(type==3) sTemp="pepxml";
	sTemp=Ext;
	lvi.iSubItem	= 2;
	lvi.pszText		= (LPTSTR)(LPCTSTR)sTemp;
	lvi.cchTextMax	= sTemp.GetLength();
	m_RList.SetItem(&lvi);
}

void CFTDRDlg::ADRList()
{
	int count=m_RList.GetItemCount();	
	for(int i=0;i<count;i++)
	{
		LV_ITEM lvi;
		lvi.mask= LVIF_TEXT; 
		CString sTemp;
		lvi.iItem		= i;
		lvi.iSubItem	= 0;
		sTemp.Format("%d",i+1);
		lvi.pszText		= (LPTSTR)(LPCTSTR)sTemp;
		lvi.cchTextMax	= sTemp.GetLength();
		m_RList.SetItem(&lvi);
	}
}

void CFTDRDlg::MsgThreadEnd(LPARAM lParam)
{
	ThreadNum--;
	prog+=detProg;
	m_WorkPro.SetPos(prog);
	Parm->ThreadID=lParam;
	if(ThreadNum==0)
	{
		(CButton *)GetDlgItem(IDC_CAL)->EnableWindow(true);	
		WriteReport();
		Parm->FreeRepBuff();
		m_WorkPro.SetPos(100);
		MessageBox("Finshed!");	
	}
}

bool CFTDRDlg::MsgTaskFinshed(LPARAM lParam)
{
	CReport rept;
	Parm->GetRepBuff(lParam,&rept);
	Parm->ResetRepState(lParam);//this message has been processed avoid another time
	OTRACE("Finish: %s\n",rept.rFile);	
	CString sMsg;
	sMsg.Format("Finish: %s\n",rept.rFile);
	AddMsg(sMsg);
	Add_Error_Msg(&rept);
	//rept=*(CReport *)wParam;rept.count>0&&
	if(rept.mcount>0)
	{
		if(rept.IsRepValidate)
		{
			RepList.push_back(rept);
			WriteReport(rept);
		}
		//else return false;
	}	
	if(CurrentFP>=ValidNum) //tail the thread end itself, and will send back a MSG_TASK_FINISHED message
	{	
		Parm->IsEnd=true;
		Parm->ThreadID=lParam;
		SetEvent(Parm->hEventP);
		return true;
	}
	//normal task renew process
	CString st1,st3,st4;
	CString rtype,dtype;
	int type;
	st1=m_RList.GetItemText(CurrentFP,1);
	st4=m_RawList.GetItemText(CurrentFP,1);
	rtype=m_RList.GetItemText(CurrentFP,2);
	st3=m_OutPath;
	//AdjustSVMFName();
	if((rtype=="htm"||rtype=="dat"))type=R_TYPE_MASCOT;//Mascot html or dat	
	else if(rtype=="Path")type=R_TYPE_SEQUEST;//Sequest output path
	else if(rtype=="pepxml") type=R_TYPE_PEPXML;//PepXML
	else if(rtype=="massmatrix") type=R_TYPE_MASSMATRIX;//massmatrix
	//no others, the add rotunie has checked the format, no poosible for others
	WaitForSingleObject(Parm->hEventE,INFINITE);
	OTRACE("Sending %d: %s!\n", CurrentFP,st1);
	sMsg.Format("Sending %d: %s!\n", CurrentFP,(LPCTSTR)st1);
	AddMsg(sMsg);
    CurrentFP++;		
	Parm->GL_CriticalSection.Lock();			
	Parm->RFile=st1;
	Parm->RawFile=st4;
	Parm->sOutPath=st3;		
	Parm->type=type;
	Parm->ThreadID=lParam;
	Parm->IsDeal=false;
	SetEvent(Parm->hEventP);
	ResetEvent(Parm->hEventE);
	Parm->GL_CriticalSection.Unlock();	
	prog+=detProg;
	m_WorkPro.SetPos(prog);
	return true;
}

//processed some noticed message, to update the process bar more preciely
void CFTDRDlg::MsgRawNotices(WPARAM wParam,LPARAM lParam)
{
	CString sMsg;
	if(wParam==MSG_RAW_LOAD)
	{
		sMsg="The raw file:";
		sMsg+=Parm->RawFile;
		sMsg+=" Load into memory finished!";
		AddMsg(sMsg);
	}
	else if(wParam==MSG_RESULT_LOAD)
	{
		sMsg="The results in file: ";
		sMsg+=Parm->RFile;
		sMsg+=" Load finished!";
		AddMsg(sMsg);
	}
	else if(wParam==MSG_RAW_LOAD)
	{
		sMsg="raw file: ";
		sMsg+=Parm->RawFile;
		sMsg+=" Load finished!";
		AddMsg(sMsg);
	}
	else if(wParam==MSG_MODEL_END)
	{
		sMsg="Building calibration model for ";
		sMsg+=Parm->RFile;
		sMsg+=" finished!";
		AddMsg(sMsg);
	}
	else if(wParam==MSG_NOIDRESULTS)
	{
		sMsg="No identification results in ";
		sMsg+=Parm->RFile;		
		AddMsg(sMsg);
	}
}

LRESULT CFTDRDlg::OnDisplay(WPARAM wParam,LPARAM lParam)
{
	int bReturn=0;
	if(wParam==MSG_THREAD_END) MsgThreadEnd(lParam);
	else if(wParam==MSG_TASK_FINISHED) bReturn=MsgTaskFinshed(lParam);
	else MsgRawNotices(wParam,lParam);	
	return bReturn;
}

BOOL CFTDRDlg::DestroyWindow() 
{
	// TODO: Add your specialized code here and/or call the base class
	delete Parm;	
	return CDialog::DestroyWindow();
}

void CFTDRDlg::WriteReport()
{
	CStdioFile dFile;
	CString sTemp;
	sTemp=m_OutPath+"\\"+"FTDR_report.txt";
	if(!dFile.Open(sTemp,CFile::modeCreate|CFile::modeWrite|CFile::typeText))
	{
		MessageBox("Can not create the report file!");
		return;
	}
	size_t count=RepList.size();
	for(size_t i=0;i<count;i++)
	{		
		dFile.WriteString("Database search file:");
		dFile.WriteString(RepList[i].rFile);
		dFile.WriteString("\n");

		sTemp.Format("The validate identifications:%d\n",RepList[i].mcount);
		dFile.WriteString(sTemp);

		dFile.WriteString("\n");
		sTemp.Format("Recalibrations:%d, and ppb level is(by MS count model):%d\n",RepList[i].count,RepList[i].ppbLevel);
		dFile.WriteString(sTemp);

		size_t k,num;

		if(Parm->ModelType==MODEL_SVM)
		{
			sTemp.Format("too many svm model parameters, not output in this file.\n");
			dFile.WriteString(sTemp);
		}
		else if(Parm->ModelType==MODEL_LIN)
		{
			sTemp.Format("The linear model parameters is:\n");
			dFile.WriteString(sTemp);
			num=RepList[i].LCoff.size();
			for(k=0;k<num;k++)
			{
				sTemp.Format(" %lf\n",RepList[i].LCoff[k]);
				dFile.WriteString(sTemp);
			}
		}

		num=RepList[i].SelParName.size();

		sTemp.Format("The %d selected features for modeling are:\n",num);
		dFile.WriteString(sTemp);
		for(k=0;k<num;k++)
		{			
			dFile.WriteString(RepList[i].SelParName[k].c_str());
			dFile.WriteString("\n");
		}

		sTemp.Format("Before calibratation: \nmean=%f\nstd=%f\nMET=%fppm\n\n",RepList[i].mtA.mean,RepList[i].mtA.std,3*RepList[i].mtA.std);
		dFile.WriteString(sTemp);

		sTemp.Format("After calibratation \nmean=%f\nstd=%f\nMET=%fppm\n\n",RepList[i].mtB.mean,RepList[i].mtB.std,3*RepList[i].mtB.std);
		dFile.WriteString(sTemp);

		sTemp.Format("MET complex model building:%d\n",RepList[i].IsMETModel);
		dFile.WriteString(sTemp);
		if(RepList[i].IsMETModel)
		{
			sTemp.Format("The parameters of complex MET model are: a=%lf\tb=%lf\tLikeHood=%lf\n\n",
				RepList[i].METModel[1],RepList[i].METModel[2],RepList[i].METModel[0]);
			dFile.WriteString(sTemp);
		}

		if(MS1OFormat.GetCurSel()!=0)
		{
			sTemp.Format("The SVM model based MS1 calibration Parameters: mean=%lf\tstd=%lf\tMET=%lf\n",
				RepList[i].Model_mean,RepList[i].Model_std,3*RepList[i].Model_std);
			dFile.WriteString(sTemp);
		}

		Write_Error_code(&RepList[i],&dFile);
	}
	dFile.Close();
}

void CFTDRDlg::WriteReport(int i)
{
	CStdioFile dFile;
	CString sTemp;
	size_t count=RepList.size();
	if(i>=count||i<=0) return;
	sTemp.Format("FTDR_report_%d.txt",i);
	sTemp=m_OutPath+"\\"+sTemp;
	if(!dFile.Open(sTemp,CFile::modeCreate|CFile::modeWrite|CFile::typeText))
	{
		MessageBox("Can not create the report file!");
		return;
	}	

	dFile.WriteString("Database search file:");
	dFile.WriteString(RepList[i].rFile);
	dFile.WriteString("\n");

	sTemp.Format("The validate identifications:%d\n",RepList[i].mcount);
	dFile.WriteString(sTemp);

	dFile.WriteString("\n");
	sTemp.Format("Recalibrations:%d, and ppb level is(by MS count model):%d\n",RepList[i].count,RepList[i].ppbLevel);
	dFile.WriteString(sTemp);	

	size_t k,num;
	if(Parm->ModelType==MODEL_SVM)
	{
		sTemp.Format("too many svm model parameters, not output in this file.\n");
		dFile.WriteString(sTemp);
	}
	else if(Parm->ModelType==MODEL_LIN)
	{
		sTemp.Format("The linear model parameters is:\n");
		dFile.WriteString(sTemp);
		num=RepList[i].LCoff.size();
		for(k=0;k<num;k++)
		{
			sTemp.Format(" %lf\n",RepList[i].LCoff[k]);
			dFile.WriteString(sTemp);
		}
	}
	
	num=RepList[i].SelParName.size();

	sTemp.Format("The %d selected features for modeling are:\n",num);
	dFile.WriteString(sTemp);
	for(k=0;k<num;k++)
	{			
		dFile.WriteString(RepList[i].SelParName[k].c_str());
		dFile.WriteString("\n");
	}

	sTemp.Format("Before calibratation: \nmean=%f\nstd=%f\nMET=%fppm\n\n",RepList[i].mtA.mean,RepList[i].mtA.std,3*RepList[i].mtA.std);
	dFile.WriteString(sTemp);

	sTemp.Format("After calibratation \nmean=%f\nstd=%f\nMET=%fppm\n\n",RepList[i].mtB.mean,RepList[i].mtB.std,3*RepList[i].mtB.std);
	dFile.WriteString(sTemp);

	sTemp.Format("MET complex model building:%d\n",RepList[i].IsMETModel);
	dFile.WriteString(sTemp);
	if(RepList[i].IsMETModel)
	{
		sTemp.Format("The parameters of complex MET model are: a=%lf\tb=%lf\tLikeHood=%lf\n\n",
			RepList[i].METModel[1],RepList[i].METModel[2],RepList[i].METModel[0]);
		dFile.WriteString(sTemp);
	}

	if(MS1OFormat.GetCurSel()!=0)
	{
		sTemp.Format("The SVM model based MS1 calibration Parameters: mean=%lf\tstd=%lf\tMET=%lf\n",
			RepList[i].Model_mean,RepList[i].Model_std,3*RepList[i].Model_std);
		dFile.WriteString(sTemp);
	}

	Write_Error_code(&RepList[i],&dFile);	
	dFile.Close();
}

void CFTDRDlg::WriteReport(CReport &rt)
{
	CStdioFile dFile;
	CString sTemp=GetPureFName(rt.RawFile);
	sTemp+="_RepFile.txt";
	sTemp=m_OutPath+"\\"+sTemp;
	if(!dFile.Open(sTemp,CFile::modeCreate|CFile::modeWrite|CFile::typeText))
	{
		MessageBox("Can not create the report file!");
		return;
	}	

	dFile.WriteString("Database search file:");
	dFile.WriteString(rt.rFile);
	dFile.WriteString("\n");

	sTemp.Format("The validate identifications:%d\n",rt.mcount);
	dFile.WriteString(sTemp);

	dFile.WriteString("\n");
	sTemp.Format("Recalibrations:%d, and ppb level is(by MS count model):%d\n",rt.count,rt.ppbLevel);
	dFile.WriteString(sTemp);	

	size_t k,num;

	if(Parm->ModelType==MODEL_SVM)
	{
		sTemp.Format("too many svm model parameters, not output in this file.\n");
		dFile.WriteString(sTemp);
	}
	else if(Parm->ModelType==MODEL_LIN)
	{
		sTemp.Format("The linear model parameters is:\n");
		dFile.WriteString(sTemp);
		num=rt.LCoff.size();
		for(k=0;k<num;k++)
		{
			sTemp.Format(" %lf\n",rt.LCoff[k]);
			dFile.WriteString(sTemp);
		}
	}

	num=rt.SelParName.size();

	sTemp.Format("The %d selected features for modeling are:\n",num);
	dFile.WriteString(sTemp);
	for(k=0;k<num;k++)
	{			
		dFile.WriteString(rt.SelParName[k].c_str());
		dFile.WriteString("\n");
	}

	sTemp.Format("Before calibratation: \nmean=%f\nstd=%f\nMET=%fppm\n\n",rt.mtA.mean,rt.mtA.std,3*rt.mtA.std);
	dFile.WriteString(sTemp);

	sTemp.Format("After calibratation \nmean=%f\nstd=%f\nMET=%fppm\n\n",rt.mtB.mean,rt.mtB.std,3*rt.mtB.std);
	dFile.WriteString(sTemp);

	sTemp.Format("MET complex model building:%d\n",rt.IsMETModel);
	dFile.WriteString(sTemp);
	if(rt.IsMETModel)
	{
		sTemp.Format("The parameters of complex MET model are: a=%lf\tb=%lf\tLikeHood=%lf\n\n",
			rt.METModel[1],rt.METModel[2],rt.METModel[0]);
		dFile.WriteString(sTemp);
	}

	if(MS1OFormat.GetCurSel()!=0)
	{
		sTemp.Format("The SVM model based MS1 calibration Parameters: mean=%lf\tstd=%lf\tMET=%lf\n",
			rt.Model_mean,rt.Model_std,3*rt.Model_std);
		dFile.WriteString(sTemp);
	}

	Write_Error_code(&rt,&dFile);	
	dFile.Close();
}

bool CFTDRDlg::Write_Error_code(CReport *rt,CStdioFile *dFile)
{
	if(dFile==NULL||rt==NULL) return false;
	if(rt->Error_code==0x00) return false;
	dFile->WriteString("Error Messages in this task:\n");
	if(rt->Error_code&ERROR_RAWFAILUR)
	{
		dFile->WriteString(ErrorMessage[0]);
		dFile->WriteString("\n");
	}
	if(rt->Error_code&ERROR_NOTLOAD)
	{
		dFile->WriteString(ErrorMessage[1]);
		dFile->WriteString("\n");
	}
	if(rt->Error_code&ERROR_COUPUTF)
	{
		dFile->WriteString(ErrorMessage[2]);
		dFile->WriteString("\n");
	}
	if(rt->Error_code&ERROR_MDFAILUR)
	{
		dFile->WriteString(ErrorMessage[3]);
		dFile->WriteString("\n");
	}
	if(rt->Error_code&ERROR_RNOPEN)
	{
		dFile->WriteString(ErrorMessage[4]);
		dFile->WriteString("\n");
	}
	if(rt->Error_code&ERROR_NOS_RFORMAT)
	{
		dFile->WriteString(ErrorMessage[5]);
		dFile->WriteString("\n");
	}
	if(rt->Error_code&ERROR_MS1OUTF)
	{
		dFile->WriteString(ErrorMessage[6]);
		dFile->WriteString("\n");
	}
	if(rt->Error_code&ERROR_MODEL)
	{
		dFile->WriteString(ErrorMessage[7]);
		dFile->WriteString("\n");
	}
	return true;
}

void CFTDRDlg::AddMsg(CString msg)
{
	//if(MsgLine>100)
	//{
	//	m_MsgList.ResetContent();
	//	MsgLine=0;
	//}
	/*m_MsgList.AddString(msg);*/
	MSGPostBox.SetSel(-1,-1);
	CString tmp=msg;
	tmp+="\r\n";
	MSGPostBox.ReplaceSel((LPCTSTR)tmp);
	//MsgLine++;
}

void CFTDRDlg::Add_Error_Msg(CReport *rt)
{
	if(rt==NULL) return;
	if(rt->Error_code==0x00) return;
	CString str="Error Messages in this task:\n";
	AddMsg(str);	
	if(rt->Error_code&ERROR_RAWFAILUR)
	{
		AddMsg(ErrorMessage[0]);
		
	}
	if(rt->Error_code&ERROR_NOTLOAD)
	{
		AddMsg(ErrorMessage[1]);		
	}
	if(rt->Error_code&ERROR_COUPUTF)
	{
		AddMsg(ErrorMessage[2]);		
	}
	if(rt->Error_code&ERROR_MDFAILUR)
	{
		AddMsg(ErrorMessage[3]);		
	}
	if(rt->Error_code&ERROR_RNOPEN)
	{
		AddMsg(ErrorMessage[4]);		
	}
	if(rt->Error_code&ERROR_NOS_RFORMAT)
	{
		AddMsg(ErrorMessage[5]);	
	}
	if(rt->Error_code&ERROR_MS1OUTF)
	{
		AddMsg(ErrorMessage[6]);
	}
	if(rt->Error_code&ERROR_MODEL)
	{
		AddMsg(ErrorMessage[7]);
	}
}

void CFTDRDlg::OnBtaddsraw() 
{
	// TODO: Add your control notification handler code here
	TCHAR szPathBuffer[_MAX_PATH];
	TCHAR szDrive[_MAX_DRIVE];
	TCHAR szDir[_MAX_DIR];

	CString sFileFilter(_T("data Files (*.raw) | *.raw ||"));
	_tsplitpath(m_sRawFile, szDrive, szDir, NULL, NULL);
	_tmakepath(szPathBuffer, szDrive, szDir, NULL, NULL);

	CFileDialog dlg(TRUE, 
					_T("raw"), 
					_T("*.raw"), 
					OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_EXTENSIONDIFFERENT | OFN_EXPLORER,
					sFileFilter,
					this );
	dlg.m_ofn.lpstrTitle = _T("Open raw File");
	dlg.m_ofn.lpstrDefExt = _T("raw");
	dlg.m_ofn.lpstrInitialDir = szPathBuffer;
    
	if(dlg.DoModal()==IDOK)
	{
         m_sRawFile=dlg.GetPathName();
		 AddRawItem(m_sRawFile);
	}	
}

void CFTDRDlg::OnBtbatchraw() 
{
	// TODO: Add your control notification handler code here
   CXSBrowseFolder browse;
   browse.SetTitle("Open the directory with *.raw files");
   char temppath[MAX_PATH];
   browse.ModifyStyle(BIF_EDITBOX|BIF_RETURNONLYFSDIRS);  
   if(browse.Show(this->GetSafeHwnd(),temppath)==2)
   {
	  HANDLE hFind;
	  WIN32_FIND_DATA fd;
	  CString strFileSpec =temppath;
	  if (strFileSpec.Right (1) != "\\") strFileSpec += "\\";
	  strFileSpec += "*.raw";
	  if ((hFind =::FindFirstFile ((LPCTSTR) strFileSpec, &fd))==INVALID_HANDLE_VALUE) 
	  { 
		  return ;
	  }
	  CWaitCursor wait;
	  do 
	  {
		  CString strFileName = (LPCTSTR) &fd.cFileName;
		  CString temp=temppath;
		  temp+="\\";
		  temp=temp+GetNTS(fd.cFileName);
		  temp.Replace(":\\\\",":\\");
		  AddRawItem(temp);		
	  } while (::FindNextFile (hFind, &fd));
	  ::FindClose (hFind);	
   } 	
}

void CFTDRDlg::OnBtremsraw() 
{
	// TODO: Add your control notification handler code here
	CString *pList;
	int count=m_RawList.GetItemCount();
	if(count<=0) return;
	pList=new CString[count];
	int i;
	for(i=0;i<count;i++)
	{		
		pList[i]=m_RawList.GetItemText(i,1);	
	}
	POSITION pos = m_RawList.GetFirstSelectedItemPosition();
	if (pos == NULL) return;		
	else
	{
		while (pos)
		{
			int nItem = m_RawList.GetNextSelectedItem(pos);
			pList[nItem].Empty();		
		}
	}
	m_RawList.DeleteAllItems();
	for(i=0;i<count;i++)
	{
		if(pList[i].IsEmpty()) continue;	
		 AddRawItem(pList[i]);
	}
	delete []pList;
	//ADDList();	
}

void CFTDRDlg::OnBtremaraw() 
{
	// TODO: Add your control notification handler code here
	m_RawList.DeleteAllItems();
	
}

void CFTDRDlg::AddRawItem(CString raw)
{
	int i=m_RawList.GetItemCount();	
	LV_ITEM lvi;
	lvi.mask= LVIF_TEXT; 
	CString sTemp;
	lvi.iItem		= i;
	lvi.iSubItem	= 0;
	sTemp.Format("%d",i+1);
	lvi.pszText		= (LPTSTR)(LPCTSTR)sTemp;
	lvi.cchTextMax	= sTemp.GetLength();
	m_RawList.InsertItem(&lvi); 

	lvi.iSubItem	= 1;
	sTemp=raw;
	lvi.pszText		= (LPTSTR)(LPCTSTR)sTemp;
	lvi.cchTextMax	= sTemp.GetLength(); 
	m_RawList.SetItem(&lvi);
}

void CFTDRDlg::OnBnClickedBtviewhist()
{
	// TODO: 在此添加控件通知处理程序代码
	HistView dlg;
	size_t count=RepList.size();
	if(count<=0) 
	{
		dlg.DoModal();
	}
	else
	{	
		HistList Ht[2];		
		for(size_t i=0;i<count;i++)
		{		
			if(RepList[i].IsRepValidate)
			{
				Ht[0].push_back(RepList[i].mtA);	
				Ht[1].push_back(RepList[i].mtB);	
			}
		}		
		dlg.SetData(Ht);
		dlg.DoModal();
	}
}

//void CFTDRDlg::OnBnClickedBtModels()
//{
//	// TODO: 在此添加控件通知处理程序代码
//	TCHAR szPathBuffer[_MAX_PATH];
//	TCHAR szDrive[_MAX_DRIVE];
//	TCHAR szDir[_MAX_DIR];
//
//	CString sFileFilter(_T("svm model files (*.txt) | *.txt ||"));
//	_tsplitpath(m_sRawFile, szDrive, szDir, NULL, NULL);
//	_tmakepath(szPathBuffer, szDrive, szDir, NULL, NULL);
//
//	CFileDialog dlg(FALSE, 
//					_T("text"), 
//					_T("*.txt"), 
//					OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_EXTENSIONDIFFERENT | OFN_EXPLORER,
//					sFileFilter,
//					this );
//	dlg.m_ofn.lpstrTitle = _T("Save svm model File");
//	dlg.m_ofn.lpstrDefExt = _T("txt");
//	dlg.m_ofn.lpstrInitialDir = szPathBuffer;
//    
//	if(dlg.DoModal()==IDOK)
//	{
//         svm_file=dlg.GetPathName();		 
//	}	
//}

void CFTDRDlg::OnBnClickedBtaddsetting()
{
	// TODO: 在此添加控件通知处理程序代码
	AdvSetDlg dlg;
	AdvPars tmpAdvpars=APars;
	dlg.SetParStore(&tmpAdvpars);
	if(dlg.DoModal()==IDOK)
	{
		if(dlg.IsValidate)
		{
			APars=tmpAdvpars;
			UpdateEffective();
		}
		else MessageBox("some parameters error!");
	}
}

void CFTDRDlg::UpdateEffective()
{
	MINRT=APars.MinRT;//普通XIC搜索中，RT范围最小值，默认-5.0min
    MAXRT=APars.MaxRT;//普通XIC搜索中，RT范围最大值，默认5.0min
	MIN_SN=APars.Min_SN;
	GD_CUT=APars.Iso_GD;//同位素峰的你和优度门限
	MAX_CH=APars.Max_Ch;//最大的考虑电荷
	INT_TIME_MAX=APars.Xic_Int_Max;//XIC构建中，信号缺失的最大次数
	RINT_CUT=APars.Min_Rel_Int;//相对信号强度过滤门限
	ISO_CUT=APars.Min_ISO;//means no cut here同位素峰数目的限制
	MAX_INT_RT=APars.RT_Max_Int;
	CH_RT_MIN=APars.Ch_RT_Min;//电荷确定中，搜索图谱叠加的区间向左扩展范围
	CH_RT_MAX=APars.Ch_RT_Max;//电荷确定中，搜索图谱叠加的区间向右扩展范围
	
	MAX_SVM_TRAIN=APars.TrainSize;//SVM模型训练样本数量最大值
	MIN_H_ISO_NUM=APars.MinIsoNum;//扩展查找母离子时，同位素峰最小值
	R_INT_CUT_H=APars.RIntCTH;//扩展查找时，相对信号强度最小值
	GD_CUT_H= APars.GdCtH;//original is 0.2
	
	LLEXTEND=APars.LExt;// 10
	RLEXTEND=APars.RExt;// 10

	PRE_DM=APars.Pre_dm;//10ppm
	XICCALIBRATE=APars.XICCalibrate;//true
	PARIONRED=APars.PIonResearch;//true;
	SVMC=APars.SVM_C;

	IS_DYA_XIC=APars.IsDyaXIC;
	DYA_XIC_SW=APars.XIC_SW;
	DYA_XIC_SM=APars.XIC_SMC;
}

void CFTDRDlg::OnBnClickedBtloadds()
{
	// TODO: 在此添加控件通知处理程序代码
	// TODO: 在此添加控件通知处理程序代码
	TCHAR szPathBuffer[_MAX_PATH];
	TCHAR szDrive[_MAX_DRIVE];
	TCHAR szDir[_MAX_DIR];

	CString sFileFilter(_T("experiment design file (*.txt) | *.txt ||"));
	_tsplitpath(m_sRawFile, szDrive, szDir, NULL, NULL);
	_tmakepath(szPathBuffer, szDrive, szDir, NULL, NULL);

	CFileDialog dlg(TRUE, 
					_T("text"), 
					_T("*.txt"), 
					OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_EXTENSIONDIFFERENT | OFN_EXPLORER,
					sFileFilter,
					this );
	dlg.m_ofn.lpstrTitle = _T("Load experiment design to File");
	dlg.m_ofn.lpstrDefExt = _T("txt");
	dlg.m_ofn.lpstrInitialDir = szPathBuffer;
    
	if(dlg.DoModal()==IDOK)
	{
         CString SFName=dlg.GetPathName();	
		 if(!LoadEXPFromFile(SFName)) MessageBox("Can not Load, check the file path and name, or file format, or parameters range!");
		 else 
		 {
			 UpdateEffective();	
			 UpdateData(FALSE);
		 }
	}	
}

void CFTDRDlg::OnBnClickedBtsaveds()
{
	// TODO: 在此添加控件通知处理程序代码
	TCHAR szPathBuffer[_MAX_PATH];
	TCHAR szDrive[_MAX_DRIVE];
	TCHAR szDir[_MAX_DIR];

	CString sFileFilter(_T("experiment design file (*.txt) | *.txt ||"));
	_tsplitpath(m_sRawFile, szDrive, szDir, NULL, NULL);
	_tmakepath(szPathBuffer, szDrive, szDir, NULL, NULL);

	CFileDialog dlg(FALSE, 
					_T("text"), 
					_T("*.txt"), 
					OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_EXTENSIONDIFFERENT | OFN_EXPLORER,
					sFileFilter,
					this );
	dlg.m_ofn.lpstrTitle = _T("Save experiment design to File");
	dlg.m_ofn.lpstrDefExt = _T("txt");
	dlg.m_ofn.lpstrInitialDir = szPathBuffer;
    
	if(dlg.DoModal()==IDOK)
	{
         CString SFName=dlg.GetPathName();	
		 if(!SaveEXPtoFile(SFName)) MessageBox("Can not save, check the file path and name!");
	}	
}


bool CFTDRDlg::SaveEXPtoFile(CString &SFName)
{
	CStdioFile dFile;
	if(!dFile.Open(SFName,CFile::modeCreate|CFile::typeText|CFile::modeWrite)) return false;
	CString tmp;
	//file paths	
	dFile.WriteString("[RAW_PATH]\n");
	int count=m_RawList.GetItemCount();
	int i;
	for(i=0;i<count;i++)
	{		
		tmp=m_RawList.GetItemText(i,1);	
		dFile.WriteString(tmp);
		dFile.WriteString("\n");
	}

	dFile.WriteString("[DBSR_PATH]\n");
	count=m_RList.GetItemCount();
	for(i=0;i<count;i++)
	{		
		tmp=m_RList.GetItemText(i,1);	
		tmp+="\t";
		tmp+=m_RList.GetItemText(i,2);
		dFile.WriteString(tmp);
		dFile.WriteString("\n");
	}

	dFile.WriteString("[OUTPATH]\n");
	dFile.WriteString(m_OutPath);
	dFile.WriteString("\n");

	dFile.WriteString("[BASE_PARS]\n");
	tmp.Format("ThreadNum=%d\n",m_ThreadNum.GetCurSel());
	dFile.WriteString(tmp);
	
	tmp.Format("Is_Out_Cal_Data=%d\n",m_OPD.GetCheck());
	dFile.WriteString(tmp);	

	tmp.Format("Is_Out_SVM=%d\n",m_RPO.GetCheck());
	dFile.WriteString(tmp);

	tmp.Format("Is_use_FT_Status_Parameters=%d",StatusPar.GetCheck());	
	dFile.WriteString(tmp);
	dFile.WriteString("\n");

	tmp.Format("Model_Type=%d\n",ModelSel.GetCurSel());
	dFile.WriteString(tmp);

	tmp.Format("MS1_OUT=%d\n",MS1OFormat.GetCurSel());
	dFile.WriteString(tmp);

	tmp.Format("MS2_OUT=%d\n",MS2Type.GetCurSel());
	dFile.WriteString(tmp);

	dFile.WriteString("[Advance_Pars]\n");
	APars.SaveToFile(dFile);
	dFile.WriteString("[Mascot_Filter]\n");
	Parm->dFilter.WriteToFile(dFile);
	dFile.WriteString("[Sequest_Filter]\n");
	Parm->dFilter1.WriteToFile(dFile);
	return true;
}

bool CFTDRDlg::LoadEXPFromFile(CString &SFName)
{
	CStdioFile dFile;
	if(!dFile.Open(SFName,CFile::typeText|CFile::modeRead)) return false;
	CString tmp;
	//file paths	
	int idx;
	CString item;
	AdvPars tmpAdvpars;

	m_RawList.DeleteAllItems();
	m_RList.DeleteAllItems();

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	if(tmp!="[RAW_PATH]") goto FAILURE_Local;

	while(dFile.ReadString(tmp))
	{
		if(tmp=="[DBSR_PATH]") break;
		AddRawItem(tmp);
	}
	
	while(dFile.ReadString(tmp))
	{
		if(tmp=="[OUTPATH]") break;
		idx=tmp.Find("\t");
		if(idx==-1) continue;
		item=tmp.Mid(0,idx);
		idx++;
		tmp=tmp.Mid(idx);
		AddRitem(item,tmp);		
	}

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	m_OutPath=tmp;

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	if(tmp!="[BASE_PARS]") goto FAILURE_Local;

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	int RD=sscanf((LPCTSTR)tmp,"ThreadNum=%d",&idx);
	if(RD!=1)goto FAILURE_Local;
	m_ThreadNum.SetCurSel(idx);


	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	RD=sscanf((LPCTSTR)tmp,"Is_Out_Cal_Data=%d",&idx);
	if(RD!=1)goto FAILURE_Local;	
	m_OPD.SetCheck(idx);

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	RD=sscanf((LPCTSTR)tmp,"Is_Out_SVM=%d",&idx);
	if(RD!=1) goto FAILURE_Local;
	m_RPO.SetCheck(idx);

	//if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	//idx=tmp.Find("SVM_Model_Path=");
	//if(idx==-1)	goto FAILURE_Local;
	//svm_file=tmp.Mid(15);
	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	RD=sscanf((LPCTSTR)tmp,"Is_use_FT_Status_Parameters=%d",&idx);
	if(RD!=1) goto FAILURE_Local;
	StatusPar.SetCheck(idx);

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	RD=sscanf((LPCTSTR)tmp,"Model_Type=%d",&idx);
	if(RD!=1)	goto FAILURE_Local;
	ModelSel.SetCurSel(idx);	

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	RD=sscanf((LPCTSTR)tmp,"MS1_OUT=%d",&idx);
	if(RD!=1) goto FAILURE_Local;
	MS1OFormat.SetCurSel(idx);

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	RD=sscanf((LPCTSTR)tmp,"MS2_OUT=%d",&idx);
	if(RD!=1) goto FAILURE_Local;
	MS2Type.SetCurSel(idx);	

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	if(tmp!="[Advance_Pars]")  goto FAILURE_Local;	
	bool BR=tmpAdvpars.LoadFromFile(dFile);
	if(tmpAdvpars.CheckValidate()) APars=tmpAdvpars;	
	else BR=false;

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	if(tmp!="[Mascot_Filter]")  goto FAILURE_Local;
	if(!Parm->dFilter.LoadFromFile(dFile)) goto FAILURE_Local;
	if(!Parm->dFilter.IsValidate())goto FAILURE_Local;

	if(!dFile.ReadString(tmp)) goto FAILURE_Local;
	if(tmp!="[Sequest_Filter]")  goto FAILURE_Local;
	if(!Parm->dFilter1.LoadFromFile(dFile)) goto FAILURE_Local;
	if(!Parm->dFilter1.IsValidate())goto FAILURE_Local;

	dFile.Close();
	return BR;
FAILURE_Local:
	dFile.Close();
	return false;
}
void CFTDRDlg::OnCbnSelchangeCombmodel()
{
	// TODO: 在此添加控件通知处理程序代码

}
