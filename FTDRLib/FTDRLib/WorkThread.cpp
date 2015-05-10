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
	while(1)//�߳�һ���������᲻�ϵȴ���������ֱ�����̷߳���MSG_THREAD_END����
	{	
		//�����Ѿ���ɣ������߳��Զ���ֹ	
		if(parm->IsEnd)
		{	
			QuantObj.PostFalseMsg(MSG_THREAD_END);
			AfxEndThread(0);			
			return 0;
		}	
		//�ȴ�����
		WaitForSingleObject(parm->hEventP,INFINITE);
		//��������ǲ����Ѿ��������������쳣�ַ������ԣ����µȴ�
		if(parm->IsDeal) continue;
		//�������ID�ǲ��Ƿ����Լ���
		//if(parm->ThreadID!=QuantObj.ThreadID) continue;
		//��ʼ������������������������ı�
		parm->GL_CriticalSection.Lock();
		parm->IsDeal=true;
		TParm=*parm;
		SetEvent(parm->hEventE);	
		ResetEvent(parm->hEventP);
		parm->GL_CriticalSection.Unlock();	
		//��ʼУ��������parm���Խ����µĲ������
		QuantObj.ExecuteCal(&TParm);
		//QuantObj.OutPutCalData(&TParm);
		//��������
		QuantObj.GenReport(&tRep);		
		tRep.rFile=TParm.RFile;	
		tRep.RawFile=TParm.RawFile;
		tRep.IsRepValidate=true;//�ñ�����Ч�����߳̿��Դӹ�������ȡ
		if(!QuantObj.IsOutSuccess) tRep.Error_code|=ERROR_MS1OUTF;
		parm->SetRepBuff(QuantObj.ThreadID,&tRep);

		printf("Msg %d: %s!\n", QuantObj.ThreadID,TParm.RFile);

		//����������ɵ���Ϣ
		QuantObj.PostFalseMsg(MSG_TASK_FINISHED);
		//T_Report(&dFile,&tRep);		
	}    
}

