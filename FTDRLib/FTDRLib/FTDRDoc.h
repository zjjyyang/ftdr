// FTDRDoc.h : interface of the CFTDRDoc class
//


#pragma once


class CFTDRDoc : public CDocument
{
protected: // create from serialization only
	CFTDRDoc();
	DECLARE_DYNCREATE(CFTDRDoc)

// Attributes
public:

// Operations
public:

// Overrides
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CFTDRDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
};


