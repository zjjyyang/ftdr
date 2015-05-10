// HistView.cpp : 实现文件
//

#include "stdafx.h"
#include "FTDR.h"
#include "HistView.h"
#include "math.h"

// HistView 对话框

IMPLEMENT_DYNAMIC(HistView, CDialog)

HistView::HistView(CWnd* pParent /*=NULL*/)
	: CDialog(HistView::IDD, pParent)
	, m_Current(0)
	, m_total(0)
	, m_rawfile(_T(""))
{

}

HistView::~HistView()
{
}

void HistView::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_ETRAW, m_rawfile);
	DDX_Text(pDX, IDC_EDIT2, m_Current);
	DDX_Text(pDX, IDC_EDIT1, m_total);
}


BEGIN_MESSAGE_MAP(HistView, CDialog)
	ON_WM_PAINT()
	ON_BN_CLICKED(IDC_BTLAST, &HistView::OnBnClickedBtlast)
	ON_BN_CLICKED(IDC_BTNEXT, &HistView::OnBnClickedBtnext)
	ON_BN_CLICKED(IDC_BTLOAD, &HistView::OnBnClickedBtload)
	ON_BN_CLICKED(IDC_BTSAVE, &HistView::OnBnClickedBtsave)
END_MESSAGE_MAP()


// HistView 消息处理程序

void HistView::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	// TODO: 在此处添加消息处理程序代码
	// 不为绘图消息调用 CDialog::OnPaint()
	CRect rect; 
	GetClientRect(rect); 
	CRect LF,RT;
	LF.top=rect.top;
	RT.top=rect.top;
	LF.bottom=rect.bottom-120;
	RT.bottom=LF.bottom;
	int W=(rect.Width()-10)/2;
	LF.right=W;
	RT.right=rect.right;
	LF.left=rect.left;
	RT.left=W+10;
	myHistogram mt;
	if(m_Current<m_total)
	{
		HistsA.at(m_Current,mt);
		DrawHist(&dc,LF,mt);
		HistsB.at(m_Current,mt);
		DrawHist(&dc,RT,mt,true);
	}
	//CDialog::OnPaint(); 
}

bool HistView::DrawHist(CPaintDC *pDC,CRect &rect,myHistogram &Hist,bool Iscenter)
{
	//CPaintDC pDC(this); // device context for painting

	CPen  *oldpen=NULL;
	CBrush *oldbrush=NULL;
	CBrush Whitebrush(RGB(255,255,255));
	CBrush Bluebrush(RGB(0,0,255));
	CPen BluePen(PS_SOLID,1,RGB(0,0,255));	
	pDC->SetTextAlign(TA_CENTER);
	CFont *oldfont,newfont;
		VERIFY(newfont.CreateFont(
		12,                        // nHeight
		0,                         // nWidth
		0,                         // nEscapement
		0,                         // nOrientation
		FW_NORMAL,                 // nWeight
		FALSE,                     // bItalic
		FALSE,                     // bUnderline
		0,                         // cStrikeOut
		ANSI_CHARSET,              // nCharSet
		OUT_DEFAULT_PRECIS,        // nOutPrecision
		CLIP_DEFAULT_PRECIS,       // nClipPrecision
		DEFAULT_QUALITY,           // nQuality
		DEFAULT_PITCH | FF_SWISS,  // nPitchAndFamily
		"Arial"));                 // lpszFacename
	oldfont=(CFont *)pDC->SelectObject(&newfont);
/////////////////////////////////////////////////////绘制坐标系
	oldbrush = (CBrush *)pDC->SelectObject(&Whitebrush);	
	pDC->FillRect(&rect,&Whitebrush);
	pDC->SelectObject(oldbrush);
	double maxy=Hist.GetMaxY();
	if(maxy<=0) return false;
	double sxMax=Hist.GetMaxX();
	double sxMin=Hist.GetMinX();
	if(Iscenter)
	{
		double stmp=fabs(sxMin);
		if(stmp>sxMax) sxMax=stmp;
		else if(sxMax>0) sxMin=-sxMax;		
	}
	if(sxMax<=sxMin) return false;

	float xscale=(rect.Width()-40.0f)/(sxMax-sxMin);
	float yscale=(rect.Height()-30)/maxy;
	pDC->MoveTo(rect.left+30,rect.top+10);
	pDC->LineTo(rect.left+30,rect.bottom-20);////////////////y axe
	int basey=rect.bottom-20;  
	pDC->MoveTo(rect.left+30,basey);
	pDC->LineTo(rect.right-20,basey);/////////////////x axe   
	int basex=rect.left+30;
	int tempy;
	CString tick;
	float bindata=maxy/10.0f;
	int i;
	for(i=0;i<11;i++)
	{
        tempy=int(basey-i*bindata*yscale);
        pDC->MoveTo(basex,tempy);
	    pDC->LineTo(rect.right-20,tempy);	          
	    tick.Format("%.4f",i*bindata);
		pDC->TextOut(basex-15,tempy-5,tick);
	}///////////////////////////////////////////////////////y axe tick
	int tempx;
    bindata=(sxMax-sxMin)/20.0f;
    float tempMAX=20;
    if(bindata<=4) 
	{
		bindata=1;
		tempMAX=sxMax-sxMin;
	}
	for(i=0;i<tempMAX;i++)
	{
		tempx=(int)(basex+i*bindata*xscale);
			//tempx=(int)(basex+(rect.Width()-40.0f)*i/20);
		pDC->MoveTo(tempx,basey-2);
	    pDC->LineTo(tempx,basey+2);	          
	    tick.Format("%4.2f",i*bindata+sxMin);
	    pDC->TextOut(tempx,basey+5,tick);
	}///////////////////////////////////////////////////////////x axe tick
	oldpen=(CPen *)pDC->SelectObject(&BluePen);		
	CPoint point1,point2;
	oldbrush = (CBrush *)pDC->SelectObject(&Bluebrush);	
	for(i=0;i<Hist.n;i++)
	{
		point1.x=(long)((Hist.range[i]-sxMin)*xscale+basex);
		point1.y=(long)(basey-Hist.bin[i]*yscale);		
		//pDC->MoveTo(point);
		point2.y=basey;	
		point2.x=(long)((Hist.range[i+1]-sxMin)*xscale+basex);
		if(point2.x<=point1.x) point2.x=point1.x+1;
		//pDC->LineTo(point);	
		pDC->Rectangle(point1.x,point1.y,point2.x,point2.y);
		CRect trect(point1,point2);
		pDC->FillRect(&trect,&Bluebrush);

	} 
	pDC->SelectObject(oldpen);
    pDC->SelectObject(oldfont);
	pDC->SelectObject(oldbrush);
	//draw the guassian PDF curve
	CPen RedPen(PS_SOLID,1,RGB(255,0,0));	
	oldpen=(CPen *)pDC->SelectObject(&RedPen);
	double bins=(sxMax-sxMin)/200;
	if(bins*xscale<1) bins=xscale;
	double t=sxMin+bins;
	point1.x=(long)((t-sxMin)*xscale+basex);
	point1.y=(long)(basey-Hist.GetProbGS(t)*yscale);	
	pDC->MoveTo(point1);
	while(t<sxMax)
	{
		t+=bins;
		point1.x=(long)((t-sxMin)*xscale+basex);
		point1.y=(long)(basey-Hist.GetProbGS(t)*yscale);
		pDC->LineTo(point1);
	}
	pDC->SelectObject(oldpen);
	return true;
}
void HistView::OnBnClickedBtlast()
{
	// TODO: 在此添加控件通知处理程序代码
	if(m_Current>0)
	{
		m_Current--;
		myHistogram mt;
		HistsA.at(m_Current,mt);
		m_rawfile=mt.rawname;
		UpdateData();
		Invalidate();
	}
}

void HistView::OnBnClickedBtnext()
{
	// TODO: 在此添加控件通知处理程序代码
	if(m_Current<m_total-1)
	{
		m_Current++;
		myHistogram mt;
		HistsA.at(m_Current,mt);
		m_rawfile=mt.rawname;
		UpdateData();
		Invalidate();
	}
}

void HistView::RemTail(char *buf)
{
	int slen=strlen(buf);
	if(slen<1) return;
	if(slen>=2&&buf[slen-2]=='\r')
	{
		buf[slen-2]='\0';
	}
	else if(buf[slen-1]=='\n')
	{
		buf[slen-1]='\0';
	}
}

bool HistView::SaveData(char *fname)
{
	FILE *fp;
	fp=fopen(fname,"w");
	if(fp==NULL) return false;
	myHistogram mt;
	char sRawName[MAX_PATH];
	int CT=HistsA.size();
	int k,i=HistsB.size();
	CT=CT>i?i:CT;
	for(i=0;i<CT;i++)
	{
		HistsA.at(i,mt);	
		scopy(sRawName,mt.rawname);
		fprintf(fp,">H0\t%d\t%lf\t%lf\t%s\n",mt.n,mt.mean,mt.std,sRawName);
		for(k=0;k<mt.n;k++)
		{
			fprintf(fp,"%lf\t%lf\n",mt.range[k],mt.bin[k]);
		}
		fprintf(fp,"%lf\t0.0\n",mt.range[mt.n]);
		HistsB.at(i,mt);
		scopy(sRawName,mt.rawname);
		fprintf(fp,">H1\t%d\t%lf\t%lf\t%s\n",mt.n,mt.mean,mt.std,sRawName);
		for(k=0;k<mt.n;k++)
		{
			fprintf(fp,"%lf\t%lf\n",mt.range[k],mt.bin[k]);
		}
		fprintf(fp,"%lf\t0.0\n",mt.range[mt.n]);
	}
	fclose(fp);
	return true;
}
bool HistView::LoadData(char *fname)
{
	FILE *fp;
	fp=fopen(fname,"r");
	if(fp==NULL) return false;
	myHistogram mt;
	char buf[1024];
	HistsB.clear();
	HistsA.clear();
	char *pStr;
	int type=-1;
	int cidx;
	char sRawName[MAX_PATH];
	fgets(buf,1023,fp);
	RemTail(buf);
	pStr=strstr(buf,">H");	
	if(pStr!=buf)
	{
		MessageBox("The histogram file format is not correct!");
		fclose(fp);
		return false;
	}
	//can not deal the finlename with space,sucha as 
	//E:\sample data\9.mzXML
	sscanf(buf,">H%d\t%d\t%lf\t%lf\t%s",&type,&mt.n,&mt.mean,&mt.std,sRawName);	
	if(mt.n>0)mt.alloc(mt.n);	
	mt.rawname=sRawName;
	cidx=0;
	while(!feof(fp))
	{
		fgets(buf,1023,fp);
		RemTail(buf);
		pStr=strstr(buf,">H");
		if(pStr==buf)
		{
			if(type==0)HistsA.push_back(mt);
			else HistsB.push_back(mt);
			sscanf(buf,">H%d\t%d\t%lf\t%lf\t%s",&type,&mt.n,&mt.mean,&mt.std,sRawName);
			if(mt.n>0)mt.alloc(mt.n);
			mt.rawname=sRawName;
			cidx=0;
			continue;
		}
		else
		{
			if(mt.n>0)
			{
				double dt1,dt2;
				sscanf(buf,"%lf\t%lf",&dt1,&dt2);
				if(cidx<mt.n)
				{
					mt.range[cidx]=dt1;
					mt.bin[cidx]=dt2;
					cidx++;
				}
				else if(cidx==mt.n)
				{
					mt.range[cidx]=dt1;
				}
			}
		}
	}
	if(type==0)HistsA.push_back(mt);
	else HistsB.push_back(mt);
	fclose(fp);	
	m_total=HistsA.size();
	int i=HistsB.size();
	m_total=m_total>i?i:m_total;
	m_Current=0;
	if(m_total>0)
	{
		HistsA.at(0,mt);	
		m_rawfile=mt.rawname;
	}
	UpdateData(false);
	Invalidate();	
	return true;
}
void HistView::OnBnClickedBtload()
{
	// TODO: 在此添加控件通知处理程序代码
	TCHAR szPathBuffer[_MAX_PATH];
	TCHAR szDrive[_MAX_DRIVE];
	TCHAR szDir[_MAX_DIR];
	CString sDataFile;

	CString sFileFilter(_T("Histogram data Files (*.hist) | *.hist ||"));
	_tsplitpath(sDataFile, szDrive, szDir, NULL, NULL);
	_tmakepath(szPathBuffer, szDrive, szDir, NULL, NULL);

	CFileDialog dlg(TRUE, 
					_T("hist"), 
					_T("*.hist"), 
					OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_EXTENSIONDIFFERENT | OFN_EXPLORER,
					sFileFilter,
					this );
	dlg.m_ofn.lpstrTitle = _T("Open hist File");
	dlg.m_ofn.lpstrDefExt = _T("hist");
	dlg.m_ofn.lpstrInitialDir = szPathBuffer;
    
	if(dlg.DoModal()==IDOK)
	{
         sDataFile=dlg.GetPathName();
		 char fname[MAX_PATH];
		 scopy(fname,sDataFile);	
		 LoadData(fname);
	}
}

void HistView::scopy(char *buf,CString str)
{
	int slen=str.GetLength();
	for(int i=0;i<slen;i++)
	{
		buf[i]=str.GetAt(i);
	}
	buf[slen]='\0';
}
void HistView::OnBnClickedBtsave()
{
	// TODO: 在此添加控件通知处理程序代码
	if(HistsA.size()<=0||HistsB.size()<=0) return;
	TCHAR szPathBuffer[_MAX_PATH];
	TCHAR szDrive[_MAX_DRIVE];
	TCHAR szDir[_MAX_DIR];
	CString sDataFile;

	CString sFileFilter(_T("Histogram data Files (*.hist) | *.hist ||"));
	_tsplitpath(sDataFile, szDrive, szDir, NULL, NULL);
	_tmakepath(szPathBuffer, szDrive, szDir, NULL, NULL);

	CFileDialog dlg(FALSE, 
					_T("hist"), 
					_T("*.hist"), 
					OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_EXTENSIONDIFFERENT | OFN_EXPLORER,
					sFileFilter,
					this );
	dlg.m_ofn.lpstrTitle = _T("Open hist File");
	dlg.m_ofn.lpstrDefExt = _T("hist");
	dlg.m_ofn.lpstrInitialDir = szPathBuffer;
    
	if(dlg.DoModal()==IDOK)
	{
         sDataFile=dlg.GetPathName();
		 char fname[MAX_PATH];
		 scopy(fname,sDataFile);	
		 SaveData(fname);
	}
}

void HistView::SetData(HistList A[2])
{
	HistsA=A[0];
	HistsB=A[1];	
	m_total=HistsA.size();
	int i=HistsB.size();
	m_total=m_total>i?i:m_total;
	m_Current=0;
	myHistogram mt;
	if(m_total>0)
	{
		HistsA.at(0,mt);
		m_rawfile=mt.rawname;
	}
	//UpdateData();
	//Invalidate();
}
