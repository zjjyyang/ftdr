// FTDRDoc.cpp : implementation of the CFTDRDoc class
//

#include "stdafx.h"
#include "FTDR.h"

#include "FTDRDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CFTDRDoc

IMPLEMENT_DYNCREATE(CFTDRDoc, CDocument)

BEGIN_MESSAGE_MAP(CFTDRDoc, CDocument)
END_MESSAGE_MAP()


// CFTDRDoc construction/destruction

CFTDRDoc::CFTDRDoc()
{
	// TODO: add one-time construction code here

}

CFTDRDoc::~CFTDRDoc()
{
}

BOOL CFTDRDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}




// CFTDRDoc serialization

void CFTDRDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}


// CFTDRDoc diagnostics

#ifdef _DEBUG
void CFTDRDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CFTDRDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CFTDRDoc commands
