/*
 *  
 *
 *  Created by Maria Abou Chakra on 24/02/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 * CR1
 */

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<math.h>
#include<map>
#include <string.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<vector>
#include<algorithm>

#include "randnumgenCR1-D.h"
#include "individualCR1-D.h"

using namespace std;
typedef map<string,int>Mymap;

int main(int argc, char* argv[]);
void outputFunction(std::vector<Individual> &popN, Mymap & phenolist,Mymap & genolist,double *totalC, double *genC,int const& roundN,int const& countGen,int const& inD, int const& oneinD,string optN, int const& gent);
double probFN(double const& prob, double t, double x,double const& gamma,double const& pF);
bool memberQ(int *rl, int const& num, int const& listlen);
void goop(int *grp, int const& gInd, int const& inD, RandNum Rnum);
void oFn(std::vector<Individual> &popN, int const& roundN, int const& inD, int const& gInd, double const& prob,RandNum RNUm, double factor,double const& inP,double const& gamma,double const& pF,double *totalC);
void gFn(std::vector<Individual> &popN, int const& roundN, int const& inD, int const& gInd, double const& sCo,int const& gameNum,double const& prob,RandNum RNUm, double factor,double const& inP,double const& gamma,double const& pF,double *totalC);
void mgFn(std::vector<Individual> &popN,int const& muprob,double *mutation, int const& inD,int const& roundN,double const& sigma,RandNum RNUm);
void wfFn(std::vector<Individual> &popN,  std::vector<Individual> &newP, int const& NinD, RandNum RNUm);

int main(int argc, char* argv[])
{
	RandNum	RNUm;
	string nameF (argv[1]);
	string loca(argv[2]);
	int runnum=atoi(argv[3]);
	string outOpt(argv[4]);
	int gameNum = atoi(argv[5]);
	int inD=atoi(argv[6]);
	double prob = atof(argv[7]);
	int muprob = atoi(argv[8]);
	double sCo=atof(argv[9]);
	int gInd = atoi(argv[10]);
	int roundN = atoi(argv[11]);
	int mPy=atoi(argv[12]);
	double targetInfo=atof(argv[13]);
	double factor=atof(argv[14]);
	int indTypeNum = atoi(argv[15]);
	int repSim=atoi(argv[16]);
	int maxgen=atoi(argv[17]);
	double inP=atof(argv[18]);
	double gamma=atof(argv[19]);
	double pF=atof(argv[20]);
	double sigma=atof(argv[21]);
	//int group[gInd];
	int i,id,num,numf,numl,fl,temp;
	int fixed = 0;
	int countGen = 0;
	double *genC = new double[4];
	Mymap phenolist;
	Mymap genolist;
	map<string,int>::iterator iter;
	string phenoN("");
	char pC[10];
	
	if(outOpt.find("4",0)!= string::npos)
	{
		
		i=0;
		for(i=0;i<(roundN); i++)
		{
			
			if(i>=(roundN/2))//someof the behaviours need to happen only half of the game
			{ 
				
				if((indTypeNum==1)){sprintf(pC, "%d", mPy);}
				else if((indTypeNum==2)){sprintf(pC, "%d", 0);}
				else if((indTypeNum==3)){sprintf(pC, "%d", mPy/2);}
				else if((indTypeNum==4)){sprintf(pC, "%d", mPy);}
				else if((indTypeNum==5)){sprintf(pC, "%d", 0);}
			}
			else{
			 	if((indTypeNum==1)){sprintf(pC, "%d", 0);}
				else if((indTypeNum==2)){sprintf(pC, "%d", mPy);}
				else if((indTypeNum==3)){sprintf(pC, "%d", mPy/2);}
				else if((indTypeNum==4)){sprintf(pC, "%d", mPy);}
				else if((indTypeNum==5)){sprintf(pC, "%d", 0);}
			}
						
		    phenoN=phenoN+pC;
			
			if(i < (roundN-1))//doesnt add a dash after the last item
			{
				phenoN=phenoN+"-";
			}
		}
	
	}
	
	ofstream fout;
	double mutation[muprob];
	i=0;
	for(i=0;i<muprob;i++)
	{
		mutation[i]=atof(argv[22+i]);
	}
	
	double * totalC = new double[((roundN)+10)];//caclulates the sum of all the columns over all generations and individuals
    //initializes the population
	std::vector<Individual> popN(inD);  
	std::vector<Individual> newP(inD);
	
	int repeatSim = 0;
	for(repeatSim=0; repeatSim<=repSim; repeatSim++)
	{
		//re initialise variables
		fixed = 0;
		countGen = 0;
		i=0;
		for(i=0;i<((roundN)+10);i++)
		{
			totalC[i]=0;
		}
		iter=phenolist.begin();
		for(iter=phenolist.begin(); iter !=phenolist.end();iter++)
		{
			phenolist[iter->first]=0;
		}
		iter=genolist.begin();
		for(iter=genolist.begin(); iter !=genolist.end();iter++)
		{
			genolist[iter->first]=0;
		}
		
		//initialize individuals and types in a population
		
		i=0;
		for(i=0;i<inD;i++)
		{	
			if(repeatSim>0){popN[i].ReInitialize(roundN);}
			else{popN[i].Initialize(roundN);}
		
			popN[i].SetValues(i, roundN, mPy,indTypeNum,RNUm);
			popN[i].SetTarget(double(gInd)*double(roundN));
			popN[i].SetIndTarget(double(gInd));
		
			
		}
				
		do
		{
			
			gFn(popN, roundN, inD, gInd, sCo,gameNum,prob,RNUm, factor, inP,gamma,pF,totalC);
			
			
			if(muprob==0)
			{
				fixed=0;
				i=0;
				for(i=1; i<inD;i++)
				{
					
					if (popN[1].GetInfo(0) != popN[i].GetInfo(0) )
					{
						fixed = 1;
						i=inD;
					}	
					
					
				}
			}
			else{fixed =1;}
			
			if((fixed<=0)||(muprob>0))
			{
				
				if(outOpt.find("4",0)!= string::npos)
				{phenolist[phenoN]=0;}
				genC[0]=0;genC[1]=0;genC[2]=0;genC[3]=0;
				i = 0; 
				for(i=0; i < inD;i++)
				{ 
					outputFunction(popN,phenolist,genolist,totalC, genC,roundN,countGen,inD, i,outOpt, maxgen);
					
				}
				
				if(outOpt.find("4",0)!= string::npos)
				{
					
					if(phenolist[phenoN]<=double(inD/2))
					{fixed=0;}
				}
				
				
				if (outOpt.find("6",0)!= string::npos)
				{
					
					nameF=argv[1];
					sprintf(pC, "%d", runnum+repeatSim);
					nameF= loca+"Gall-"+nameF+" "+pC+"-a.dat";
					char *filN6 = new char[(nameF.length()+1)];
					strcpy(filN6, nameF.c_str());
					
					fout.open(filN6,ios::app); 
					if (!fout) 
					{
						//cout << "Unable to open" << argv[1] << "for appending.\n"; 
						fout.open(filN6); // open for writing
					}
					fout<<genC[0]/double(inD);fout<<"  ";
					fout<<genC[1]/double(inD);fout<<"  ";
					fout<<genC[2]/double(inD);fout<<"  ";
					fout<<genC[3]/double(inD);fout<<"  ";
					
					fout<<"\n";
					fout.close(); 
				}
				
			}
			
			
			wfFn(popN,  newP, inD, RNUm);
			mgFn(popN,muprob,mutation, inD,roundN, sigma,RNUm);
			
			
			i = 0; 
			for(i=0; i < inD;i++)
			{
				popN[i].ResetValues();
			}
			countGen++;
			
			if((muprob>0)&&(countGen==maxgen))
			{fixed=0;}
			
		}while(fixed==1);
		
		if (outOpt.find("1",0)!= string::npos)
		{
			
			nameF=argv[1];
			nameF= loca+"tall-"+nameF+"-a.dat";
			char *filN1 = new char[(nameF.length()+1)];
			strcpy(filN1, nameF.c_str());
			
			fout.open(filN1,ios::app); 
			if (!fout) 
			{
				//cout << "Unable to open" << argv[1] << "for appending.\n"; 
				fout.open(filN1); // open for writing
			}
			id=0;
			for(id=0;id<(roundN+7);id++)
			{
		
				fout<<totalC[id]/double(double(inD)*double(maxgen))<<"  ";
			}
						
			fout<<"\n";
			
			fout.close(); 
		}
		if (outOpt.find("2",0)!= string::npos)
		{
			nameF=argv[1];
			nameF= loca+"fpc-"+nameF+"-a.dat";
			char *filN2 = new char[(nameF.length()+1)];
			strcpy(filN2, nameF.c_str());
			fout.open(filN2,ios::app); 
			if (!fout) 
			{
				//cout << "Unable to open" << argv[1] << "for appending.\n"; 
				fout.open(filN2); // open for writing
			}
			
			iter=phenolist.begin();
			for(iter=phenolist.begin(); iter !=phenolist.end();iter++)
			{
				fout<<iter->first<<" "<<iter->second<<" "<<runnum+repeatSim<<"\n";
			}
			
			fout.close(); 
		}	
		if (outOpt.find("3",0)!= string::npos)
		{
			nameF=argv[1];
			sprintf(pC, "%d", runnum+repeatSim);
			nameF= loca+"fgc-"+nameF+"-"+pC+"-a.dat";
			char *filN3 = new char[(nameF.length()+1)];
			strcpy(filN3, nameF.c_str());
			fout.open(filN3,ios::app); 
			if (!fout) 
			{
				//cout << "Unable to open" << argv[1] << "for appending.\n"; 
				fout.open(filN3); // open for writing
			}
			
			iter=genolist.begin();
			for(iter=genolist.begin(); iter !=genolist.end();iter++)
			{
				fout<<iter->first<<" "<<iter->second<<" "<<runnum+repeatSim<<"\n";
			}
			
			fout.close(); 
		}	
		if(outOpt.find("4",0)!= string::npos)
		{
			nameF=argv[1];
			nameF= loca+"stability-"+nameF+"-a.dat";
			char *filN4 = new char[(nameF.length()+1)];
			strcpy(filN4, nameF.c_str());
			fout.open(filN4,ios::app); 
			if (!fout) 
			{
				//cout << "Unable to open" << argv[1] << "for appending.\n"; 
				fout.open(filN4); // open for writing
			}
			
			fout<<(countGen)<<" "<<phenolist[phenoN]<<" "<<phenoN<<"\n";
			fout.close();
		}	
		
		
	}
	
	///memory release
	delete [] totalC;
	delete [] genC;
	genC = NULL;
	totalC = NULL;
	
	return 0;
}
void outputFunction(std::vector<Individual> &popN,Mymap  &phenolist,Mymap &genolist,double *totalC, double *genC,int const& roundN,int const& countGen,int const& inD, int const& oneinD,string optN, int const& gent)
{
	int i,id,num,numf,numl;
	double contribution;
	char numstr[200];
	string phenotype;
	string genotype;
	
	if (optN.find("1",0)!= string::npos)
	{
		if(popN[oneinD].GetInfo(2)>0){totalC[0]=totalC[0]+(double(popN[oneinD].GetInfo(3))/double(popN[oneinD].GetInfo(2)));}
		if(popN[oneinD].GetInfo(2)>0){totalC[1]=totalC[1]+(double(popN[oneinD].GetGoalReached())/double(popN[oneinD].GetInfo(2)));}
		totalC[2]=totalC[2]+popN[oneinD].GetEndowment();
	}
	id=0;
	num=4;

		id=0;
		contribution=0;
		for(id=0;id<(roundN);id++)
		{
			
			if(num==4)
			{	
				if(popN[oneinD].GetInfo(2)>0)
				{
					if((optN.find("1",0)!= string::npos)||(optN.find("6",0)!= string::npos)){contribution=contribution+double(popN[oneinD].GetB(id));}
					if((optN.find("2",0)!= string::npos)||(optN.find("4",0)!= string::npos)||(optN.find("3",0)!= string::npos))
					{
						sprintf(numstr, "%d", popN[oneinD].GetB(id)/popN[oneinD].GetInfo(2));
						phenotype=phenotype+numstr;
					}
				}
				else
				{
					if((optN.find("2",0)!= string::npos)||(optN.find("4",0)!= string::npos)||(optN.find("3",0)!= string::npos)){phenotype=phenotype+"-";}
					if(optN.find("1",0)!= string::npos){contribution=-1;}
				}
				
				if(id < (roundN-1))
				{
					if((optN.find("2",0)!= string::npos)||(optN.find("4",0)!= string::npos)||(optN.find("3",0)!= string::npos)){phenotype=phenotype+"-";}
				}
				
			}				   
			
		}	
			
		if((num==4)&&(optN.find("6",0)!= string::npos))
		{
			genC[0]=genC[0]+countGen;
			genC[1]=genC[1]+contribution/double(popN[oneinD].GetInfo(2));
			genC[2]=genC[2]+(double(popN[oneinD].GetInfo(3))/double(popN[oneinD].GetInfo(2)));
			genC[3]=genC[3]+(double(popN[oneinD].GetGoalReached())/double(popN[oneinD].GetInfo(2)));
		}
		
		if((num==4)&&(optN.find("1",0)!= string::npos))
		{
			
			if((contribution/double(popN[oneinD].GetInfo(2)))==0){totalC[roundN+3]=totalC[roundN+3]+1;}
			else if((contribution/double(popN[oneinD].GetInfo(2)))<popN[oneinD].GetIndTarget()){totalC[roundN+4]=totalC[roundN+4]+1;}
			else if((contribution/double(popN[oneinD].GetInfo(2)))==popN[oneinD].GetIndTarget()){totalC[roundN+5]=totalC[roundN+5]+1;}
			else if((contribution/double(popN[oneinD].GetInfo(2)))>popN[oneinD].GetIndTarget()){totalC[roundN+6]=totalC[roundN+6]+1;}
			
		}

	if((optN.find("2",0)!= string::npos)||(optN.find("4",0)!= string::npos))
	{phenolist[phenotype]++;}
	
	if(optN.find("3",0)!= string::npos)
	{
		sprintf(numstr, "%d ", countGen);
		genotype=numstr+phenotype;
		genolist[genotype]++;
	}
	return;
}
bool memberQ(int *rl, int const& num, int const& listlen)
{
	bool opt=false;
	int m;
	m = std::count (rl, rl+listlen, num);
	if(m>0){opt = true;}
	return opt;
}
void goop(int *grp, int const& gInd, int const& inD, RandNum Rnum)
{
	int i=0; 
	int num;
	for(i=0; i<gInd;)
	{
		num = Rnum.GetRand(0,(inD));
		grp[i]=num;
		if (i>0)
		{
			if (memberQ(grp, num, (i))==false)
			{i++;}
			
		}		
		else {i++;}
		
	}
	return;
}
double probFN(double const& prob, double t, double x,double const& gamma,double const& pF)
{
	return (prob)*(1/((pF)+exp(gamma*(x-t))));
}
void oFn(std::vector<Individual> &popN, int const& roundN, int const& inD, int const& gInd, double const& prob, RandNum RNUm, double factor,double const& inP, double const& gamma,double const& pF,double *totalC)
{
	
	int i, r;
	int num=0;
	double temp;
	double totalA=0.0;
	double totalpR=0.0;
	int totalPopIn=0;
	
	int group[gInd];
	
	goop(group, gInd, inD, RNUm);
	
	i = 0; 
	for(i=0; i < inD;i++)
	{ 
		popN[i].SetInfo(1,0);
	}
	
	r=0;
	for(r=0; r<roundN; r++)
	{
		i = 0; 
		for(i=0; i < gInd;i++)
		{	
			
			num=popN[group[i]].CalcPayment(r,totalpR);
			if (r>0)
			{
				if(num >(popN[group[i]].GetEndowment()-popN[group[i]].GetInfo(1)))
				{num=(popN[group[i]].GetEndowment()-popN[group[i]].GetInfo(1));}
			}
			popN[group[i]].SetB(r, (popN[group[i]].GetB(r)+num));
			totalC[3+r]=totalC[3+r]+double(num);
			popN[group[i]].SetInfo(1,(popN[group[i]].GetInfo(1)+num));
			totalA = totalA + num;
		}
		totalA=totalA + (totalpR * inP);
		totalpR=totalA;
	}
	
	int goalY=1;
   double	target = popN[0].GetTarget();
	
	temp=0;
	if (factor > 0)
	{
		do
		{
			temp=RNUm.GetGaussian(factor, target);
		}while (((temp<0)||(temp>(gInd*popN[0].GetEndowment()))));
		target=temp;
	}
	
	if(totalpR<target) 
	{   
		goalY=0;
	}
	num=1;
	if((pF==1)||(goalY==0))
	{
		if(gamma>=0)
		{
			
			if (RNUm.GetRand01()<=probFN(prob, target, totalpR,gamma,pF))
			{num=0;}
		}
		if(gamma<0)
		{
			num=0;
			if (RNUm.GetRand01()>probFN(prob, target, totalpR,gamma,pF))
			{num=0;}
		}
	}
	i=0;
	for(i=0; i < gInd;i++)
	{
		totalPopIn=totalPopIn + popN[group[i]].GetInfo(1);
	}
	i = 0; 
	for(i=0; i < gInd;i++)
	{
		popN[group[i]].SetInfo(2,(popN[group[i]].GetInfo(2)+1));
		popN[group[i]].SetGoalReached(popN[group[i]].GetGoalReached()+goalY);
		popN[group[i]].Payoff(num);
		
	}
	return;
}
void gFn(std::vector<Individual> &popN,  int const& roundN, int const& inD, int const& gInd, double const& sCo,int const& gameNum,double const& prob,RandNum RNUm, double factor,double const& inP,double const& gamma,double const& pF,double *totalC)
{

	for(int i=0; i<gameNum;i++)
	{
		
		oFn(popN, roundN, inD, gInd, prob,RNUm, factor, inP, gamma,pF,totalC);
	}
	
	for(int i=0; i < inD;i++)
	{
		if(popN[i].GetInfo(2)!=0)
		{
			popN[i].CalcFitness(sCo);
		}
	}
	return;
}
void mgFn(std::vector<Individual> &popN,int const& muprob,double *mutation, int const& inD,int const& roundN,double const& sigma,RandNum RNUm)
{

	for(int n=0; n<inD;n++)
	{
		for(int i=0; i<muprob; i++)
		{
			if(RNUm.GetRand01()<=mutation[i])
			{
				popN[n].SetInfo(0, popN[n].GetInfo(0)+10000);
				popN[n].Mutator(i,RNUm.GetRand(0,roundN),RNUm, sigma);
			}
		}
	}
	return;
}
double wfq(Individual player, double freqT)
{
	double wf;
	wf=((player.GetFitness()/freqT));
	return wf;
}
void wfFn(std::vector<Individual> &popN,  std::vector<Individual> &newP, int const& NinD, RandNum RNUm)
{	
	double num;
	int i,n,r,a;
	double	freqT=0;
	double wfreq[NinD];
	i=0;
	for(i=0; i<NinD; i++)
	{
		freqT = popN[i].GetFitness()+freqT;
	}
	i=0;
	wfreq[0]=wfq(popN[0], freqT);
	for(i=1; i<NinD; i++)
	{
		wfreq[i]=wfq(popN[i], freqT)+wfreq[i-1];
	}
	n=0;
	for(n=0; n<NinD;n++)
	{
		num=RNUm.GetRand01();
		i=0;
		for(i=0; i<NinD; i++)
		{
			if(num<=(wfreq[i]))
			{
				newP[n]=popN[i];
				i=NinD+1;
			}
		}
	}
	n=0;
	for(n=0; n<NinD;n++)
	{
		popN[n]=newP[n];
	}
	return;
}