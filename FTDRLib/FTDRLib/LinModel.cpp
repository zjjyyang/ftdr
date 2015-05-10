#include "StdAfx.h"
#include "LinModel.h"

LinModel::LinModel(void)
{
	COFF=NULL;
	PARNUM=0;
	stats[0]=0;
	stats[1]=0;
	stats[2]=0;
}

LinModel::~LinModel(void)
{
	if(PARNUM>0&&COFF!=NULL) delete []COFF;

}

LinModel::LinModel(const LinModel &cp)
{ 
	PARNUM=cp.PARNUM;
	if(PARNUM>0) COFF=new double[PARNUM];
	else COFF=NULL;
	for(int i=0;i<PARNUM;i++) COFF[i]=cp.COFF[i];
	stats[0]=cp.stats[0];
	stats[1]=cp.stats[1];
	stats[2]=cp.stats[2];
}

LinModel &LinModel::operator=(const LinModel &cp)
{	
	if(PARNUM>0&&COFF!=NULL) delete []COFF;
	PARNUM=cp.PARNUM;
	if(PARNUM>0) COFF=new double[PARNUM];
	else COFF=NULL;
	for(int i=0;i<PARNUM;i++) COFF[i]=cp.COFF[i];
	stats[0]=cp.stats[0];
	stats[1]=cp.stats[1];
	stats[2]=cp.stats[2];
	return *this;
}

double LinModel::Predict(double *fea,int FEATURE_NUM)
{
	double bRT=COFF[0];
	for(int i=0;i<FEATURE_NUM;i++)
	{
		bRT+=COFF[i+1]*fea[i];	
	}
	bRT=fea[0]*1e6/(1e6+bRT);//to cal back the accurate mz	
	bRT-=stats[0];
	return bRT;	
}

void LinModel::Predict(CalData &dt,vector<double> &ppmD,vector<size_t> &ParSel)
{	
	int FixedNum=ParSel.size();
	int FEATURE_NUM=FixedNum+FIXED_PAR_NUM_LIN;	
	double *Fea;
	Fea=new double[FEATURE_NUM];
	if(dt.ch==CH_UNKNOWN)
	{
		dt.PackagePar(Fea,ParSel);			
		double bRT=COFF[0];
		for(int k=0;k<FEATURE_NUM;k++)bRT+=COFF[k+1]*Fea[k];
		bRT-=stats[0];
		bRT=dt.Isomz[0]*1e6/(1e6+bRT);//to cal back the accurate mz	
		ppmD.push_back(bRT);		
	}
	else
	{
		dt.PackageParFixed(Fea,ParSel);
		for(int i=0;i<dt.IsoNum;i++)
		{
			dt.PackageParVal(Fea,i);
			double bRT=COFF[0];
			for(int k=0;k<FEATURE_NUM;k++)bRT+=COFF[k+1]*Fea[k];
			bRT-=stats[0];
			bRT=dt.Isomz[i]*1e6/(1e6+bRT);//to cal back the accurate mz
			bRT-=ISO_DIFF[i]/dt.ch;
			ppmD.push_back(bRT);	
		}
	}
	delete []Fea;
}

void LinModel::Predict(CalData &dt,vector<double> &ppmD,vector<double> &weight,vector<size_t> &ParSel)
{
	int FixedNum=ParSel.size();
	int FEATURE_NUM=FixedNum+FIXED_PAR_NUM_LIN;	
	double *Fea;
	Fea=new double[FEATURE_NUM];
	if(dt.ch==CH_UNKNOWN)
	{
		dt.PackagePar(Fea,ParSel);			
		double bRT=COFF[0];
		for(int k=0;k<FEATURE_NUM;k++)bRT+=COFF[k+1]*Fea[k];
		bRT-=stats[0];
		bRT=dt.Isomz[0]*1e6/(1e6+bRT);//to cal back the accurate mz	
		ppmD.push_back(bRT);
		weight.push_back(sqrt(dt.IsotopicE[0]));	
	}
	else
	{
		dt.PackageParFixed(Fea,ParSel);
		for(int i=0;i<dt.IsoNum;i++)
		{
			dt.PackageParVal(Fea,i);
			double bRT=COFF[0];
			for(int k=0;k<FEATURE_NUM;k++)bRT+=COFF[k+1]*Fea[k];
			bRT-=stats[0];
			bRT=dt.Isomz[i]*1e6/(1e6+bRT);//to cal back the accurate mz
			bRT-=ISO_DIFF[i]/dt.ch;
			ppmD.push_back(bRT);
			weight.push_back(sqrt(dt.IsotopicE[i]));		
		}
	}
	delete []Fea;
}


double LinModel::Predict(CalData &dt,vector<size_t> &ParSel)
{
	int FixedNum=ParSel.size();
	int FEATURE_NUM=FixedNum+FIXED_PAR_NUM_LIN;	
	double *Fea;
	Fea=new double[FEATURE_NUM];
	dt.PackagePar(Fea,ParSel);
	double bRT=COFF[0];
	for(int i=0;i<FEATURE_NUM;i++)
	{
		bRT+=COFF[i+1]*Fea[i];	
	}
	bRT-=stats[0];
	bRT=dt.Isomz[0]*1e6/(1e6+bRT);//to cal back the accurate mz
	delete []Fea;
	return bRT;	
}

int LinModel::scanSize(vector<CalData> &DataList)
{
	size_t n=DataList.size();
	int total=0;
	for(size_t i=0;i<n;i++)	total+=DataList[i].IsoNum;
	return total;
}

bool LinModel::ModelTrain(vector<CalData> &DataList,vector<size_t> &ParSel)
{
	size_t np=DataList.size();
	size_t n=scanSize(DataList);
	if(n<=3) return false;
	int FixedFea=ParSel.size();
	int FEATURE_NUM=FixedFea+FIXED_PAR_NUM_LIN;	
	//initial the parameters
	if(PARNUM>0&&COFF!=NULL) delete []COFF;
	PARNUM=FEATURE_NUM+1;
	COFF=new double[PARNUM];
	gsl_matrix *X;
	gsl_vector *y,*c;
	X=gsl_matrix_alloc(n,PARNUM);
	y=gsl_vector_alloc(n);
	c=gsl_vector_alloc(PARNUM);	
	int k,m;	
	size_t i,j=0;
	//for each item
	double *Fea;	
	Fea=new double[FEATURE_NUM];
	for(i=0;i<np;i++)
	{
		//for each iso,all the iso has the same output
		DataList[i].PackageParFixed(Fea,ParSel);
		for(k=0;k<DataList[i].IsoNum;k++)
		{
			double outputy=DataList[i].mzt+ISO_DIFF[k]/DataList[i].ch;
			outputy=(DataList[i].Isomz[k]-outputy)*1e6/outputy;
			DataList[i].PackageParVal(Fea,k);
			gsl_vector_set(y,j,outputy);
			gsl_matrix_set(X,j,0,1.0);	
			//for each fea value
			for(m=0;m<FEATURE_NUM;m++)gsl_matrix_set(X,j,m+1,Fea[m]);	
			j++;
		}
	}	
	delete []Fea;
	gsl_robustfit(X,y,c,stats);
	for(i=0;i<PARNUM;i++) COFF[i]=gsl_vector_get(c,i);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(X);	
	return true;
}