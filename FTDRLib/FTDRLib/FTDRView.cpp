// FTDRView.cpp : implementation of the CFTDRView class
//

#include "stdafx.h"
#include "FTDR.h"

#include "FTDRDoc.h"
#include "FTDRView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CFTDRView

IMPLEMENT_DYNCREATE(CFTDRView, CView)

BEGIN_MESSAGE_MAP(CFTDRView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
END_MESSAGE_MAP()

// CFTDRView construction/destruction

CFTDRView::CFTDRView()
{
	// TODO: add construction code here

}

CFTDRView::~CFTDRView()
{
}

BOOL CFTDRView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CFTDRView drawing

void CFTDRView::OnDraw(CDC* /*pDC*/)
{
	CFTDRDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: add draw code for native data here
}


// CFTDRView printing

BOOL CFTDRView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CFTDRView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CFTDRView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}


// CFTDRView diagnostics

#ifdef _DEBUG
void CFTDRView::AssertValid() const
{
	CView::AssertValid();
}

void CFTDRView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CFTDRDoc* CFTDRView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CFTDRDoc)));
	return (CFTDRDoc*)m_pDocument;
}
#endif //_DEBUG


// CFTDRView message handlers
