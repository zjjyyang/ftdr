// alQuant.cpp: implementation of the CalQuant class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "alQuant.h"
#include "OTrace.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CalQuant::CalQuant(UINT ID)
{
	ThreadID=ID;
	RWCal=NULL;
	IsOutSuccess=false;
}

CalQuant::~CalQuant()
{	
	if(RWCal!=NULL) delete RWCal;
	
}

void CalQuant::GenMF(ThreadParm *parm,CString &CalModelDataFile)
{	
	GetFNameMd(parm->RawFile,CalModelDataFile);
	CalModelDataFile=parm->sOutPath+"\\"+CalModelDataFile;
}

void CalQuant::ExecuteCal(ThreadParm *parm)
{
	if(parm->ModelType==MODEL_SVM) ExecuteCal_SVM(parm);
	else if(parm->ModelType==MODEL_LIN) ExecuteCal_Lin(parm);
	else if(parm->ModelType==MODEL_NON) ExecuteCal_Non(parm);	
}
//this function is ued for the extend calibration, added on 2010.10.18 by zhangjy
//also, you can select the XIC global calibration of simple calibration by using the parameter:parm->Model_Type;
void CalQuant::ExecuteCal_SVM(ThreadParm *parm)
{	
	//如果已经读取了一个raw文件的内容，先卸载数据
	if(RWCal!=NULL) delete RWCal;
	RWCal=new RawData;	
	//RWCal->MET=PRE_DM;
	if(IS_DYA_XIC) RWCal->InitialDFTModel();
	RWCal->IsUseExtendPar=parm->IsUseStatsPar;
	RWCal->ModelType=MODEL_SVM;
	//初始化raw文件读取接口，并且load MS1数据和MS2头
	if(!RWCal->Initial(parm->RawFile.GetBuffer()))
	{
		delete RWCal;
		parm->RawFile.ReleaseBuffer();		
		return;
	}
	parm->RawFile.ReleaseBuffer();		
	//初始化搜库结果过滤参数
	RWCal->InitialFilter(&parm->dFilter,&parm->dFilter1);
	vector<Pre_calData> IDData;
	//发送消息，通知主线程，raw文件load完毕
    PostFalseMsg(MSG_RAW_LOAD);
	//初始化校正建模数据记录文件（如果需要）
	if(parm->IsOpCData) 
	{
		CString CalModelDataFile;
		GenMF(parm,CalModelDataFile);//added on 2010.11.3
		RWCal->EnableOutput(CalModelDataFile.GetBuffer());
		CalModelDataFile.ReleaseBuffer();
	}
	else 
	{
		RWCal->IsOutput=false;
		RWCal->OIntR_Data=false;
	}

	//load搜库结果，按照过滤标准
	if(parm->type==R_TYPE_MASCOT)RWCal->CollectDataM(parm->RFile,IDData);
	else if(parm->type==R_TYPE_SEQUEST)RWCal->CollectDataS(parm->RFile,IDData);
	else if(parm->type==R_TYPE_PEPXML)	RWCal->CollectDataP(parm->RFile,IDData);
	else if(parm->type==R_TYPE_MASSMATRIX) RWCal->CollectDataX(parm->RFile,IDData);
	else RWCal->SetErrorCode(ERROR_NOS_RFORMAT);
	//通知主线程，搜库结果load完毕
	PostFalseMsg(MSG_RESULT_LOAD);
	if(IDData.size()<=0)
	{
		RWCal->SetErrorCode(ERROR_NOTLOAD);
		PostFalseMsg(MSG_NOIDRESULTS);
	}
	else
	{
		//确定输出信号强度相关模型的数据
		CString sTemp;
		GetFName(parm->RawFile,sTemp);
		if(parm->IsOpCData)
		{				
			sTemp=parm->sOutPath+"\\"+sTemp+"_IntRData.txt";	
			RWCal->SetOutputRData(sTemp.GetBuffer());
			sTemp.ReleaseBuffer();
		}
		//开始建立SVM校正模型
		if(RWCal->ModelBuilding_dMET(IDData))
		{
			//通知主线程，建模完毕
			PostFalseMsg(MSG_MODEL_END);
			RWCal->ms2format=parm->ms2format;

			//生成校正后数据存储的mgf文件名				
			GetFName(parm->RawFile,sTemp);
			sTemp=parm->sOutPath+"\\"+sTemp+".mgf";				
			//开始校正，会同步记录需要记录的数据
			RWCal->Calibrate_SVM(sTemp.GetBuffer());			
			sTemp.ReleaseBuffer();
			if(parm->IsOPSVMModel)
			{
				int slen=sTemp.GetLength();
				sTemp.SetAt(slen-1,'m');
				sTemp.SetAt(slen-2,'v');
				sTemp.SetAt(slen-3,'s');				
				RWCal->OutputModel(sTemp.GetBuffer());
				sTemp.ReleaseBuffer();
			}
			
			//add the MS1 output to mzXML function
			string tmps=parm->sOutPath.GetBuffer();
			parm->sOutPath.ReleaseBuffer();
			if(parm->ms1format!=MS1_NON)	IsOutSuccess=RWCal->OutPutMS1(tmps);
			//the center of the SVM model has been removed			
			RWCal->mtB.mean=0;			
		}
	}
	//PostFalseMsg(MSG_TASK_FINISHED);
}

//this function is ued for the extend calibration, added on 2010.10.18 by zhangjy
//also, you can select the XIC global calibration of simple calibration by using the parameter:parm->Model_Type;
void CalQuant::ExecuteCal_Lin(ThreadParm *parm)
{	
	//如果已经读取了一个raw文件的内容，先卸载数据
	if(RWCal!=NULL) delete RWCal;
	RWCal=new RawData;	
	//RWCal->MET=15;
	if(IS_DYA_XIC) RWCal->InitialDFTModel();
	RWCal->IsUseExtendPar=parm->IsUseStatsPar;	
	RWCal->ModelType=MODEL_LIN;
	//初始化raw文件读取接口，并且load MS1数据和MS2头
	if(!RWCal->Initial(parm->RawFile.GetBuffer()))
	{
		delete RWCal;
		parm->RawFile.ReleaseBuffer();		
		return;
	}
	parm->RawFile.ReleaseBuffer();		
	//初始化搜库结果过滤参数
	RWCal->InitialFilter(&parm->dFilter,&parm->dFilter1);
	vector<Pre_calData> IDData;
	//发送消息，通知主线程，raw文件load完毕
    PostFalseMsg(MSG_RAW_LOAD);
	//初始化校正建模数据记录文件（如果需要）
	if(parm->IsOpCData) 
	{
		CString CalModelDataFile;
		GenMF(parm,CalModelDataFile);//added on 2010.11.3
		RWCal->EnableOutput(CalModelDataFile.GetBuffer());
		CalModelDataFile.ReleaseBuffer();
	}
	else 
	{
		RWCal->IsOutput=false;
		RWCal->OIntR_Data=false;
	}

	//load搜库结果，按照过滤标准
	if(parm->type==R_TYPE_MASCOT)RWCal->CollectDataM(parm->RFile,IDData);
	else if(parm->type==R_TYPE_SEQUEST)RWCal->CollectDataS(parm->RFile,IDData);
	else if(parm->type==R_TYPE_PEPXML)	RWCal->CollectDataP(parm->RFile,IDData);
	else if(parm->type==R_TYPE_MASSMATRIX) RWCal->CollectDataX(parm->RFile,IDData);
	else RWCal->SetErrorCode(ERROR_NOS_RFORMAT);
	//通知主线程，搜库结果load完毕
	PostFalseMsg(MSG_RESULT_LOAD);
	if(IDData.size()<=0)
	{
		RWCal->SetErrorCode(ERROR_NOTLOAD);
		PostFalseMsg(MSG_NOIDRESULTS);
	}
	else
	{
		//确定输出信号强度相关模型的数据
		CString sTemp;
		GetFName(parm->RawFile,sTemp);
		if(parm->IsOpCData)
		{				
			sTemp=parm->sOutPath+"\\"+sTemp+"_IntRData.txt";	
			RWCal->SetOutputRData(sTemp.GetBuffer());
			sTemp.ReleaseBuffer();
		}
		//开始建立SVM校正模型
		if(RWCal->ModelBuilding_Lin(IDData))
		{
			//通知主线程，建模完毕
			PostFalseMsg(MSG_MODEL_END);
			RWCal->ms2format=parm->ms2format;
			//生成校正后数据存储的mgf文件名				
			GetFName(parm->RawFile,sTemp);
			sTemp=parm->sOutPath+"\\"+sTemp+".mgf";				
			//开始校正，会同步记录需要记录的数据
			RWCal->Calibrate_Lin(sTemp.GetBuffer());			
			sTemp.ReleaseBuffer();				
			//add the MS1 output to mzXML function
			string tmps=parm->sOutPath.GetBuffer();
			parm->sOutPath.ReleaseBuffer();
			if(parm->ms1format!=MS1_NON) IsOutSuccess=RWCal->OutPutMS1(tmps);
			//the center of the SVM model has been removed			
			RWCal->mtB.mean=0;			
		}
	}
	//PostFalseMsg(MSG_TASK_FINISHED);
}


//this function is ued for the extend calibration, added on 2010.10.18 by zhangjy
//also, you can select the XIC global calibration of simple calibration by using the parameter:parm->Model_Type;
void CalQuant::ExecuteCal_Non(ThreadParm *parm)
{	
	//如果已经读取了一个raw文件的内容，先卸载数据
	if(RWCal!=NULL) delete RWCal;
	RWCal=new RawData;
	//RWCal->MET=PRE_DM;
	if(IS_DYA_XIC) RWCal->InitialDFTModel();
	RWCal->ModelType=MODEL_NON;
	CString sTemp;
	//初始化raw文件读取接口，并且load MS1数据和MS2头
	if(!RWCal->Initial(parm->RawFile.GetBuffer()))
	{
		delete RWCal;
		parm->RawFile.ReleaseBuffer();		
		return;
	}
	parm->RawFile.ReleaseBuffer();	
	//发送消息，通知主线程，raw文件load完毕
    PostFalseMsg(MSG_RAW_LOAD);
		//初始化搜库结果过滤参数
	RWCal->InitialFilter(&parm->dFilter,&parm->dFilter1);
	vector<Pre_calData> IDData;
	//load搜库结果，按照过滤标准
	if(parm->type==R_TYPE_MASCOT)RWCal->CollectDataM(parm->RFile,IDData);
	else if(parm->type==R_TYPE_SEQUEST)RWCal->CollectDataS(parm->RFile,IDData);
	else if(parm->type==R_TYPE_PEPXML)	RWCal->CollectDataP(parm->RFile,IDData);
	else if(parm->type==R_TYPE_MASSMATRIX) RWCal->CollectDataX(parm->RFile,IDData);
	else RWCal->SetErrorCode(ERROR_NOS_RFORMAT);
	//通知主线程，搜库结果load完毕
	PostFalseMsg(MSG_RESULT_LOAD);
	if(IDData.size()<=0)
	{
		RWCal->SetErrorCode(ERROR_NOTLOAD);
		PostFalseMsg(MSG_NOIDRESULTS);
	}
	else
	{
		//确定输出信号强度相关模型的数据
		CString sTemp;
		GetFName(parm->RawFile,sTemp);
		if(parm->IsOpCData)
		{				
			sTemp=parm->sOutPath+"\\"+sTemp+"_IntRData.txt";	
			RWCal->SetOutputRData(sTemp.GetBuffer());
			sTemp.ReleaseBuffer();
		}
		//开始模型
		if(RWCal->ModelBuilding_Non(IDData))
		{
			PostFalseMsg(MSG_MODEL_END);
			RWCal->ms2format=parm->ms2format;
			//生成校正后数据存储的mgf文件名	
			GetFName(parm->RawFile,sTemp);
			sTemp=parm->sOutPath+"\\"+sTemp+".mgf";	
			//开始校正，会同步记录需要记录的数据
			RWCal->Calibrate_Non(sTemp.GetBuffer());
			sTemp.ReleaseBuffer();	
			//add the MS1 output to mzXML function
			string tmps=parm->sOutPath.GetBuffer();
			parm->sOutPath.ReleaseBuffer();
			if(parm->ms1format!=MS1_NON) IsOutSuccess=RWCal->OutPutMS1(tmps);
			//PostFalseMsg(MSG_TASK_FINISHED);
		}
	}
}

void CalQuant::OutPutCalData(ThreadParm *parm)
{	
	//如果已经读取了一个raw文件的内容，先卸载数据
	if(RWCal!=NULL) delete RWCal;
	RWCal=new RawData;	
	if(IS_DYA_XIC) RWCal->InitialDFTModel();
	//RWCal->MET=PRE_DM;
	//初始化raw文件读取接口，并且load MS1数据和MS2头
	if(!RWCal->Initial(parm->RawFile.GetBuffer()))
	{
		delete RWCal;
		parm->RawFile.ReleaseBuffer();		
		return;
	}
	parm->RawFile.ReleaseBuffer();		
	//初始化搜库结果过滤参数
	RWCal->InitialFilter(&parm->dFilter,&parm->dFilter1);
	vector<Pre_calData> IDData;
	//发送消息，通知主线程，raw文件load完毕
    PostFalseMsg(MSG_RAW_LOAD);
	//初始化校正建模数据记录文件（如果需要）
	if(parm->IsOpCData) 
	{
		CString CalModelDataFile;
		GenMF(parm,CalModelDataFile);//added on 2010.11.3
		RWCal->EnableOutput(CalModelDataFile.GetBuffer());
		CalModelDataFile.ReleaseBuffer();
	}

	//load搜库结果，按照过滤标准
	if(parm->type==R_TYPE_MASCOT)RWCal->CollectDataM(parm->RFile,IDData);
	else if(parm->type==R_TYPE_SEQUEST)RWCal->CollectDataS(parm->RFile,IDData);
	else if(parm->type==R_TYPE_PEPXML)	RWCal->CollectDataP(parm->RFile,IDData);
	else if(parm->type==R_TYPE_MASSMATRIX) RWCal->CollectDataX(parm->RFile,IDData);
	else RWCal->SetErrorCode(ERROR_NOS_RFORMAT);
	//通知主线程，搜库结果load完毕
	PostFalseMsg(MSG_RESULT_LOAD);

	//生成校正后数据存储的mgf文件名
	CString sTemp;
	GetFName(parm->RawFile,sTemp);
	sTemp=parm->sOutPath+"\\"+sTemp+"_mData.txt";	
	//开始输出SVM建模数据
	RWCal->OutPutCalData(IDData,sTemp.GetBuffer());
	sTemp.ReleaseBuffer();
	//通知主线程，建模完毕
	PostFalseMsg(MSG_MODEL_END);
}

void CalQuant::PostFalseMsg(WPARAM w_Report)
{
	WPARAM wParam=w_Report;//初始化消息参数
	LPARAM lParam=ThreadID;  
	CFrameWnd *pFrame; 
    pFrame =(CFrameWnd*)AfxGetApp()->m_pMainWnd;//Get the Current Mainframe pointer
    pFrame->PostMessage(WM_DISPLAY,wParam,lParam);//send message
}

void CalQuant::GetFName(CString fpath, CString &fname)
{
	int sidx=fpath.ReverseFind('\\');
	sidx++;	
	fname=fpath.Mid(sidx);
	sidx=fname.ReverseFind('.');
	if(sidx!=-1) fname=fname.Mid(0,sidx);
}

void CalQuant::GetFNameC(CString fpath, CString &fname)
{
	int sidx=fpath.ReverseFind('\\');	
	sidx++;
	fname=fpath.Mid(sidx);
	sidx=fname.ReverseFind('.');	
	if(sidx==-1) fname+=".cal";
	else
	{
		fname=fname.Mid(0,sidx);
		fname+=".cal";
	}
}

void CalQuant::GetFNameF(CString fpath, CString &fname)
{
	int sidx=fpath.ReverseFind('\\');	
	sidx++;
	fname=fpath.Mid(sidx);
	sidx=fname.ReverseFind('.');	
	if(sidx==-1) fname+=".fea";
	else
	{
		fname=fname.Mid(0,sidx);
		fname+=".fea";
	}
}

void CalQuant::GetFNameMd(CString fpath, CString &fname)
{
	int sidx=fpath.ReverseFind('\\');	
	sidx++;
	fname=fpath.Mid(sidx);
	sidx=fname.ReverseFind('.');	
	if(sidx==-1) fname+=".md";
	else
	{
		fname=fname.Mid(0,sidx);
		fname+=".md";
	}
}

void CalQuant::GenReport(CReport *rep)
{
	if(RWCal==NULL) rep->Error_code|=ERROR_MODEL;
	else RWCal->GenReport(rep);	
}

void CalQuant::SaveSVMModel(char *fname)
{
	RWCal->OutputModel(fname);
}

