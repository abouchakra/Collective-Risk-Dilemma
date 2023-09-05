#include <iostream>
#include<sys/time.h>
#include<cstdlib>
#include <cmath>

class RandNum
{
public:
	RandNum();//constructor fn
	RandNum(int ti);//constructor fn
	~RandNum();//destructor fn
	//int GetRandC(double min, double max); //accessor fn
	int GetRand(int min, int max); //accessor fn
	double GetRand01();
	double gasdev();//pg 289 using the Box-Muller method
	double GetGaussian( double sigma /* standard deviation */ , double mean /* mean */);
	double GetRandlog01();
	double GetRandlog(int min, int max);
private:
	int numberR;
	struct timeval tv;

};

//constructor
RandNum::RandNum()
{
	gettimeofday(&tv,NULL);
	srand(tv.tv_usec+tv.tv_sec);//seed up to the micro sec

}

//Destructor, has no action
RandNum::~RandNum()
{
}

int RandNum::GetRand(int min, int max)//just a c rand generator if needed
{
	numberR = ((double(rand())/double(RAND_MAX))*(max-min))+min;
	//numberR = (min + (max*rand()/(RAND_MAX)));
	return numberR;
}

double RandNum::GetRand01()//picks numbers betweeen 0 and 10000 but then converts it to numbers between 0 and 1
{
	return (double(GetRand(0, 1000000000))/1000000000.00);
}
double RandNum:: gasdev() /* numerical recipes C basic routine at page 289 , only this time the M twister is incorporated*/
{
	static int iset = 0;
	static double gset;
	double fac,rsq,v1,v2;
	
	if  (iset == 0) 
	{
		do {
			
			/*//v1 = 2.0* ran3() - 1.0; // picks two uniform numbers in the square from -1 to 1
			//v2 = 2.0* ran3() - 1.0;*/ //the test used ran3/ran1 but I used the M twister
			v1 = 2.0* GetRand01()- 1.0; // picks two uniform numbers in the square from -1 to 1
			v2 = 2.0* GetRand01() - 1.0;
			rsq = v1*v1 + v2*v2; //checks if they exist on the unit circle
			
		} while (rsq >= 1.0 || rsq == 0.0); //if they are not picks new ones
		
		fac = sqrt( -2.0*log( rsq )/rsq ); //Box Muller transformation
		gset = v1*fac;
		iset = 1;
		
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}
double RandNum::GetGaussian( double sigma /* standard deviation */ , double mean /* mean */)
{
	return sigma*gasdev() + mean;//in this case sigma is the interval number you wish you increase it by, it is not relative
}
