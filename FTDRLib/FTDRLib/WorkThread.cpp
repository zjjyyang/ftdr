#include "stdafx.h"
#include "WorkThread.h"
#include "alQuant.h"
//#include "OTrace.h"


UINT OutSumThread(LPVOID pParam)
{
	ThreadParm *parm;
	parm=(ThreadParm*)pParam;//Parameters
    CalQuant QuantObj(parm->ThreadID);
	SetEvent(parm->hEventB);

	printf("Create Thread %d!\n", QuantObj.ThreadID);

	CReport tRep;	
	ThreadParm TParm;
	while(1)//线程一旦创建，会不断等待接受任务，直到主线程发出MSG_THREAD_END命令
	{	
		//任务已经完成，所有线程自动终止	
		if(parm->IsEnd)
		{	
			QuantObj.PostFalseMsg(MSG_THREAD_END);
			AfxEndThread(0);			
			return 0;
		}	
		//等待任务
		WaitForSingleObject(parm->hEventP,INFINITE);
		//检查任务是不是已经处理过，如果是异常分发，忽略，重新等待
		if(parm->IsDeal) continue;
		//检查任务ID是不是发给自己的
		//if(parm->ThreadID!=QuantObj.ThreadID) continue;
		//初始化参数，共享数据区不允许改变
		parm->GL_CriticalSection.Lock();
		parm->IsDeal=true;
		TParm=*parm;
		SetEvent(parm->hEventE);	
		ResetEvent(parm->hEventP);
		parm->GL_CriticalSection.Unlock();	
		//开始校正工作，parm可以接受新的参数填充
		QuantObj.ExecuteCal(&TParm);
		//QuantObj.OutPutCalData(&TParm);
		//产生报告
		QuantObj.GenReport(&tRep);		
		tRep.rFile=TParm.RFile;	
		tRep.RawFile=TParm.RawFile;
		tRep.IsRepValidate=true;//该报告有效，主线程可以从共享区读取
		if(!QuantObj.IsOutSuccess) tRep.Error_code|=ERROR_MS1OUTF;
		parm->SetRepBuff(QuantObj.ThreadID,&tRep);

		printf("Msg %d: %s!\n", QuantObj.ThreadID,TParm.RFile);

		//发送任务完成的消息
		QuantObj.PostFalseMsg(MSG_TASK_FINISHED);
		//T_Report(&dFile,&tRep);		
	}    
}

