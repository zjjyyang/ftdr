#include "stdafx.h"
#include "SVMModel.h"
#include "OutlierRem.h"
#include "OTrace.h"
#include "math.h"
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics.h>

SVMModel::SVMModel(void)
{
	SVM_MODEL=NULL;
	svm_x_space=NULL;
	svm_prob.x=NULL;
	svm_prob.l=0;
	svm_prob.y=NULL;
	Res_mean=0;
	Res_std=0;
	gd=0;
	//C_par=64;
}

SVMModel::~SVMModel(void)
{
	if(SVM_MODEL!=NULL) svm_free_and_destroy_model(&SVM_MODEL);
	if(svm_x_space!=NULL)
	{
		free((void *)svm_prob.y);
		free((void *)svm_prob.x);
		free((void *)svm_x_space);
		svm_destroy_param(&svm_param);
	}
}

void SVMModel::Scale(int Item_Num,int FEATURE_NUM,struct svm_node *svm_x_space)
{
	size_t i,k;
	feature_max.clear();
	feature_min.clear();
	for(k=0;k<FEATURE_NUM;k++)
	{
		feature_max.push_back(svm_x_space[k].value);
		feature_min.push_back(svm_x_space[k].value);		
	}
	for(i=1;i<Item_Num;i++)
	{
		int base=i*FEATURE_NUM+i;
		for(k=0;k<FEATURE_NUM;k++)
		{
			int tmpIDX=base+k;
			if(feature_max[k]<svm_x_space[tmpIDX].value) 
				feature_max[k]=svm_x_space[tmpIDX].value;
			else if(feature_min[k]>svm_x_space[tmpIDX].value) 
				feature_min[k]=svm_x_space[tmpIDX].value;			
		}
	}

	for(k=0;k<FEATURE_NUM;k++)
	{
		if(feature_max[k]<=feature_min[k]) continue;
		double mid=feature_max[k]-feature_min[k];
		for(i=0;i<Item_Num;i++)
		{
			int tmpIDX=i*FEATURE_NUM+i+k;			
			if(svm_x_space[tmpIDX].value<=feature_min[k]) svm_x_space[tmpIDX].value=-1;
			else if(svm_x_space[tmpIDX].value>=feature_max[k])svm_x_space[tmpIDX].value=1;
			else svm_x_space[tmpIDX].value =-1+2*(svm_x_space[tmpIDX].value-feature_min[k])/mid;
		}
	}
}

//default size agreement
void SVMModel::Scale_P(svm_node *x)
{
	size_t FEATURE_NUM=feature_min.size();
	for(size_t i=0;i<FEATURE_NUM;i++)
	{
		if(x[i].value<=feature_min[i]) x[i].value=-1;
		else if(x[i].value>=feature_max[i])x[i].value=1;
		else x[i].value =-1+2*(x[i].value-feature_min[i])/\
			(feature_max[i]-feature_min[i]);
	}
}

void SVMModel::Scale_P(svm_node *x,int FEATURE_NUM)
{	
	for(size_t i=0;i<FEATURE_NUM;i++)
	{
		if(x[i].value<=feature_min[i]) x[i].value=-1;
		else if(x[i].value>=feature_max[i])x[i].value=1;
		else x[i].value =-1+2*(x[i].value-feature_min[i])/\
			(feature_max[i]-feature_min[i]);
	}
}

void SVMModel::Predict(CalData &dt,vector<double> &ppmD,vector<size_t> &ParSel)
{	
	int FixedNum=ParSel.size();
	int FEATURE_NUM=FixedNum+FIXED_PAR_NUM_SVM;
	svm_node *x;
	x=new svm_node[FEATURE_NUM+1];
	x[FEATURE_NUM].index=-1;
	int i;
	for(i=0;i<FEATURE_NUM;i++) x[i].index=i+1;	
	if(dt.ch==CH_UNKNOWN)
	{
		dt.PackagePar(x,ParSel);	
		Scale_P(x,FEATURE_NUM);
		double bRT=svm_predict(SVM_MODEL,x);
		bRT-=Res_mean;
		bRT=dt.Isomz[0]*1e6/(1e6+bRT);//to cal back the accurate mz
		ppmD.push_back(bRT);		
	}
	else
	{
		for(i=0;i<dt.IsoNum;i++)
		{		
			dt.PackagePar(x,i,ParSel);	
			Scale_P(x,FEATURE_NUM);
			double bRT=svm_predict(SVM_MODEL,x);
			bRT-=Res_mean;
			bRT=dt.Isomz[i]*1e6/(1e6+bRT);//to cal back the accurate mz
			bRT-=ISO_DIFF[i]/dt.ch;
			ppmD.push_back(bRT);
		}
	}
	delete []x;
}

void SVMModel::Predict(CalData &dt,vector<double> &ppmD,vector<double> &weight,vector<size_t> &ParSel)
{	
	int FixedNum=ParSel.size();
	int FEATURE_NUM=FixedNum+FIXED_PAR_NUM_SVM;
	svm_node *x;
	x=new svm_node[FEATURE_NUM+1];
	x[FEATURE_NUM].index=-1;
	int i;
	for(i=0;i<FEATURE_NUM;i++) x[i].index=i+1;
	if(dt.ch==CH_UNKNOWN)
	{
		dt.PackagePar(x,ParSel);	
		Scale_P(x,FEATURE_NUM);
		double bRT=svm_predict(SVM_MODEL,x);
		bRT-=Res_mean;
		bRT=dt.Isomz[0]*1e6/(1e6+bRT);//to cal back the accurate mz
		ppmD.push_back(bRT);
		weight.push_back(sqrt(dt.IsotopicE[0]));
	}
	else
	{
		for(i=0;i<dt.IsoNum;i++)
		{		
			dt.PackagePar(x,i,ParSel);	
			Scale_P(x,FEATURE_NUM);
			double bRT=svm_predict(SVM_MODEL,x);
			bRT-=Res_mean;
			bRT=dt.Isomz[i]*1e6/(1e6+bRT);//to cal back the accurate mz
			bRT-=ISO_DIFF[i]/dt.ch;
			ppmD.push_back(bRT);
			weight.push_back(sqrt(dt.IsotopicE[i]));
		}
	}
	delete []x;
}

double SVMModel::Predict(double *fea,int FEATURE_NUM)
{	
	svm_node *x;
	x=new svm_node[FEATURE_NUM+1];
	for(int i=0;i<FEATURE_NUM;i++)
	{
		x[i].index=i+1;	
		x[i].value=fea[i];
	}
	x[FEATURE_NUM].index=-1;
	Scale_P(x,FEATURE_NUM);
	double bRT=svm_predict(SVM_MODEL,x);
	bRT-=Res_mean;
	bRT=fea[0]*1e6/(1e6+bRT);//to cal back the accurate mz	
	delete []x;
	return bRT;	
}

double SVMModel::Predict(CalData &dt,vector<size_t> &ParSel)
{
	int FixedNum=ParSel.size();
	int FEATURE_NUM=FixedNum+FIXED_PAR_NUM_SVM;
	svm_node *x;
	x=new svm_node[FEATURE_NUM+1];
	x[FEATURE_NUM].index=-1;	
	for(int i=0;i<FEATURE_NUM;i++) x[i].index=i+1;
	dt.PackagePar(x,ParSel);
	Scale_P(x,FEATURE_NUM);
	double bRT=svm_predict(SVM_MODEL,x);
	bRT-=Res_mean;
	bRT=dt.Isomz[0]*1e6/(1e6+bRT);//to cal back the accurate mz
	return bRT;	
}

int SVMModel::scanSize(vector<CalData> &DataList)
{
	size_t n=DataList.size();
	int total=0;
	for(size_t i=0;i<n;i++)	total+=DataList[i].IsoNum;
	return total;
}

bool SVMModel::Model_Train(vector<CalData> &DataList,vector<size_t> &ParSel)
{	
	size_t np=DataList.size();
	if(np<=0) return false;
	int FixedFea=ParSel.size();
	int FEATURE_NUM=FixedFea+FIXED_PAR_NUM_SVM;
	int n=scanSize(DataList);
	if(n<=3) return false;
	//initial the parameters
	if(svm_x_space!=NULL)
	{
		free((void *)svm_prob.y);
		free((void *)svm_prob.x);
		free((void *)svm_x_space);
		svm_destroy_param(&svm_param);
	}
	svm_x_space=NULL;
	svm_prob.x=NULL;
	svm_prob.l=n;
	svm_prob.y=NULL;
	svm_param.svm_type = EPSILON_SVR;
	svm_param.kernel_type = RBF;
	svm_param.degree = 3;
	svm_param.gamma = 5.0/FEATURE_NUM;	// 1/k
	svm_param.coef0 = 0;
	svm_param.nu = 0.8;
	svm_param.cache_size = 1000;
	svm_param.C = SVMC;
	svm_param.eps = 1e-5;//may be modified for better performance
	svm_param.p = 0.1;
	svm_param.shrinking = 0;
	svm_param.probability = 0;
	svm_param.nr_weight = 0;
	svm_param.weight_label = NULL;
	svm_param.weight = NULL;
	//end of initial parameters
	if(SVM_MODEL!=NULL) svm_free_and_destroy_model(&SVM_MODEL);
	

	int elements=n*(FEATURE_NUM+1);
	//svm_prob.l=n;	
	svm_prob.y = (double *)malloc(sizeof(double)*n);
	svm_prob.x = (struct svm_node **)malloc(sizeof(struct svm_node *)*n);
	svm_x_space =(struct svm_node *)malloc(sizeof(struct svm_node)*elements);

	int k,m,j=0;	
	size_t i;
	int row=0;
	for(i=0;i<np;i++)
	{
		//for each iso,all the iso has the same output				
		for(k=0;k<DataList[i].IsoNum;k++)
		{
			double outputy=DataList[i].mzt+ISO_DIFF[k]/DataList[i].ch;
			outputy=(DataList[i].Isomz[k]-outputy)*1e6/outputy;	
			DataList[i].PackagePar(svm_x_space+j,k,ParSel);//a history bug for 2012.3.5			
			svm_prob.x[row] = &svm_x_space[j];
			svm_prob.y[row] = outputy;						
			for(m=0;m<FEATURE_NUM;m++)svm_x_space[j+m].index=m+1;
			j+=FEATURE_NUM;
			svm_x_space[j].index=-1;
			j++;		
			row++;
		}	
	}	

	Scale(n,FEATURE_NUM,svm_x_space);	
	SVM_MODEL = svm_train(&svm_prob,&svm_param);
	//SVM_MODEL->l=1;
	gd=CalRsquare(2);

	RMSD(n,FEATURE_NUM,svm_x_space,&svm_prob);	
	//svm_save_model(model_file_name,model);	
	return true;
}

//to estimate the error bounds
bool SVMModel::RMSD(int n,int FEATURE_NUM,struct svm_node *svm_x_space,svm_problem *svm_prob)
{
	if(n<=0) return false;
	vector<double> Res;
	for(int i=0;i<n;i++)
	{
		int tmpIDX=i*FEATURE_NUM+i;
		double fres=svm_predict(SVM_MODEL,svm_x_space+tmpIDX)-svm_prob->y[i];
		Res.push_back(fres);	
	}	
	OutlierRem ot;
	ot.ORemMixModel(Res);
	Res_mean=ot.m_new;
	Res_std=ot.std_new;
	return true;
}

//-1 means failure,0 means sucessful
bool SVMModel::OutputModel(const char *fname)
{
	if(SVM_MODEL==NULL) return false;
	return svm_save_model(fname,SVM_MODEL)+1;
}


double SVMModel::CalRsquare(int nr_fold)
{
	if(svm_prob.l<=0) return 0;
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target;
	target=new double[svm_prob.l];

	svm_cross_validation(&svm_prob,&svm_param,nr_fold,target);
	for(i=0;i<svm_prob.l;i++)
	{
		double y = svm_prob.y[i];
		double v = target[i];
		total_error += (v-y)*(v-y);
		sumv += v;
		sumy += y;
		sumvv += v*v;
		sumyy += y*y;
		sumvy += v*y;
	}
	
	//printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
	double coff=((svm_prob.l*sumvy-sumv*sumy)*(svm_prob.l*sumvy-sumv*sumy))/
			((svm_prob.l*sumvv-sumv*sumv)*(svm_prob.l*sumyy-sumy*sumy));
	printf("Cross Validation Squared correlation coefficient = %g\n",coff);
	delete []target;
	return coff;
}
