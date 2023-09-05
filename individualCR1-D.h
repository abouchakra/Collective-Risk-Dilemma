/*
 *  individualCR1.h
 *  
 *
 *  Created by Maria Abou Chakra on 2/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include<vector>
class Individual
{
public:
	Individual();//constructor fn
	void Initialize(int roundN);//constructor fn
	void ReInitialize(int roundN);//constructor fn
	~Individual();//destructor fn
	void SetValues(int num, int roundN, int maxP, int inType,RandNum RNUm);
	void ResetValues();
	int GetInfo(int r); //accessor fn
	int SetInfo(int r, int value); //accessor fn
	int GetB(int r);
	int SetB(int r,int value);
	int Getgenebelow(int r);
	int Getgeneabove(int r);
	int GetT1(int r);
	int Payoff(int &prob);
	void CalcFitness(double const& selectionCo);
	void SetFitness(double value);
	double GetFitness();
	void SetGoalReached(int value);
	int GetGoalReached();
	int GetArraySize();
	int GetEndowment();
	void SetTarget(double value);
	void SetIndTarget(double value);
	double GetTarget();
	double GetIndTarget();
	void SetMaxPayment(int value);
	int GetMaxPayment();
	void Mutator(int mp, int r,RandNum RNUm, double sigma);
	int CalcPayment(int roundN, double const & totalPotsize);
private:
	int ind[4];
	std::vector<int> bR;
	std::vector<int> genethreshold1;
	std::vector<int> genebelow;
	std::vector<int> geneabove;
	double indFit;
	int endowment;
	int maxpayment;
	double target;
	double individualtarget;
	int indTypeNum;
	int indGoalR;
	int info;
	int acol;
	int f,h;
	int flag;
};

//constructor
Individual::Individual()
{
	info=4;
	indFit=1.0;
	endowment=0;
	target=0;
	indGoalR=0;
	
}
//constructor
void Individual::Initialize(int roundN)
{
	info=4;//arraysize for individual information
	acol = roundN;
	bR.reserve(acol);
	genethreshold1.reserve(acol);
	genebelow.reserve(acol);
	geneabove.reserve(acol);
	
}
//constructor
void Individual::ReInitialize(int roundN)
{
	info=4;//arraysize for individual information
	acol = roundN;
	bR.clear();
	genethreshold1.clear();
	genebelow.clear();
	geneabove.clear();
	
}
//Destructor
Individual::~Individual()
{

}
void Individual::SetValues(int num, int roundN,int maxP, int inType,RandNum RNUm)
{
	int i, rN;
	maxpayment=maxP;
	indFit=1.0;
	target=0;
	indGoalR=0;
	indTypeNum=inType;
	endowment=2*roundN;
	ind[0]=num;
	ind[1]=0;
	ind[2]=0;
	ind[3]=0;
	i=0;
	for( i=0; i < (roundN);i++)
	{ 
		
		bR.push_back(0);//initialize the action for each round
		
		if(indTypeNum==7)
		{
			genethreshold1.push_back(RNUm.GetRand(0,1001));
			genebelow.push_back(RNUm.GetRand(0,maxpayment+1));
			geneabove.push_back(RNUm.GetRand(0,maxpayment+1));
		}
		else if(indTypeNum==6)
        {
			genethreshold1.push_back(0);
			genebelow.push_back(RNUm.GetRand(0,maxpayment+1));
            geneabove.push_back(RNUm.GetRand(0,maxpayment+1));
		}
		else if(indTypeNum==1)
		{
			genethreshold1.push_back(RNUm.GetRand(0,1001));
			if(i>=(roundN/2))
			{
				genebelow.push_back(maxpayment);
				geneabove.push_back(maxpayment);
			}
			else
			{
				genebelow.push_back(0);
				geneabove.push_back(0);
			}
		}
		else if(indTypeNum==2)
		{
			genethreshold1.push_back(RNUm.GetRand(0,1001));
			if(i>=(roundN/2))
			{
				genebelow.push_back(0);
				geneabove.push_back(0);
			}
			else
			{
				genebelow.push_back(maxpayment);
				geneabove.push_back(maxpayment);
			}
		}
		else if(indTypeNum==3)
		{
			genethreshold1.push_back(RNUm.GetRand(0,1001));
			genebelow.push_back((maxpayment/2));
			geneabove.push_back((maxpayment/2));
		}
		else if(indTypeNum==4)
		{
			genethreshold1.push_back(RNUm.GetRand(0,1001));
			genebelow.push_back((maxpayment));
			geneabove.push_back((maxpayment));
		}
		else if(indTypeNum==5)
		{
			genethreshold1.push_back(RNUm.GetRand(0,1001));
			genebelow.push_back(0);
			geneabove.push_back(0);
		}
	}	
	return;
}
void Individual::ResetValues()
{
	
	indFit=1.0;
	indGoalR=0;
	ind[1]=0;
	ind[2]=0;
	ind[3]=0;
	int i=0;
	for(i=0; i < (acol);i++)
	{ 	
		bR[i]=0;
	}
	return;
}
int Individual::GetArraySize()
{
	return info;
}
int Individual::GetInfo(int r)
{
	return ind[r];
}
int Individual::SetInfo(int r, int value)
{
	ind[r] = value;
	return ind[r];
}
int Individual::GetB(int r)
{
	return bR[r];
}
int Individual::SetB(int r,int value)
{
	bR[r] = value;
	return bR[r];
}
int Individual::Getgenebelow(int r)
{
	return genebelow[r];
}
int Individual::Getgeneabove(int r)
{
	return geneabove[r];
}
int Individual::GetT1(int r)
{
	return genethreshold1[r];
}
int Individual::Payoff(int &prob)
{
	ind[3]=ind[3]+((endowment-ind[1])*prob);
	return ind[3];
}
void Individual::CalcFitness(double const& selectionCo)
{
	
	indFit=exp(selectionCo*(double(ind[3])/double(ind[2])));
	return;
}
void Individual::SetFitness(double value)
{
	indFit=value;
	return;
}
double Individual::GetFitness()
{
	
	return indFit;
}
void Individual::SetGoalReached(int value)
{
	indGoalR = value;
	return;
}
int Individual::GetGoalReached()
{
	return indGoalR;
}
int Individual::GetEndowment()
{
	return endowment;
}
void Individual::SetTarget(double value)
{
	target=value;
	return;
}
void Individual::SetIndTarget(double value)
{
	individualtarget=double(double(target)/double(value));
	return;
}
double Individual::GetTarget()
{
	return target;
}
double Individual::GetIndTarget()
{
	return individualtarget;
}
void Individual::SetMaxPayment(int value)
{
	maxpayment = value;
	return;
}
int Individual::GetMaxPayment()
{
	return maxpayment;
}
void Individual::Mutator(int mp, int r,RandNum RNUm, double sigma)
{
	
	int temp=0;
	if(mp==0)
	{	
		genebelow[r]=RNUm.GetRand(0,maxpayment+1);
	}
	if(mp==1)
	{	
		if(indTypeNum==6)
        {
		temp=0;
		}
		else{
		do
		{temp=int(RNUm.GetGaussian(sigma, double(genethreshold1[r])/1000)*1000);
		}while (((temp<0)||(temp>1000)));
		}genethreshold1[r]=temp;
	}
	if(mp==2)
	{	
		geneabove[r]=RNUm.GetRand(0,maxpayment+1);
	}


	return;
	
}
int Individual::CalcPayment(int roundN, double const & totalPotsize)
{
	if(indTypeNum<=7)
	{
		if(((double(genethreshold1[roundN])/1000)*target)<= totalPotsize)
		{return geneabove[roundN];}
		else 
		{return genebelow[roundN];}
	}
}
