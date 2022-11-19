#pragma once
#include <vector>
#include <algorithm>
#include "Global.h"
using namespace std;

//Global variables
int Jobs[5] = { 2, 4, 6, 8, 10 };  //{20，40，60，80，100}
int Stages[4] = { 2, 3, 4 };
int Types[4] = { 1,2,3,4 };
int Jobs1[5] = { 20, 40, 60, 80,100 };
int Stages1[4] = { 3, 5, 8,10 };
int speedlevel[4] = { 1,2,3,4 };
//int SetupTime[4] = { 25, 49, 99, 124 };

//for the current schedule information
int pJob;
int pStage;
int pType;
vector<int> pLot;    //每个Job(Lot)的总批量               Jobs
vector<vector<int>> pSublot;  // 每个Lot在各阶段的分批   Jobs
vector<int> pMachines;   //每阶段的机器数                    Stages
vector<vector<double>> machineSpeed;  //the machine speed coefficient stage-speed level
vector<vector<int>> bTime;            //每个阶段基础的加工时间
vector<vector<vector<double>>> pUnitTime;  //每阶段的单位加工时间        Stages-Jobs-speed level
vector<vector<vector<int>>>pPower;//            stages-job-speed
vector<vector<int>> pSetupTime; //序列相关的启动时间 sequence-dependent setup times
vector<vector<int>> pTransferTime;  //子批在相邻阶段的转移时间 Stages-Jobs
const int MaxSublotQuantity = 5;
int DecodingType = 1;  //1: 子批优先  2： 批次优先
int MaxMS, MinMS;
int MaxTEC, MinTEC;
int ip;   //the machine power at stand-by mode
int sp;   //the machine power at setup mode

//generate instances randomly
void GenerateInstances(int InJob, int InStage, int InType, int Seed)
{
	srand(Seed);

	pJob = InJob;
	pStage = InStage;
	pType = InType;
	bTime.clear();
	pUnitTime.clear();
	pSetupTime.clear();
	pMachines.resize(pStage);
	for (int i = 0; i < pStage; i++)
	{
		pMachines[i] = 1 + rand() % 5;  //range of machines [1-5]
	}
	/*if (pType == 1)
	{
		for (int i = 0; i < pStage; i++)
		{
			if (i == 0)
			{
				pMachines[i] = 1;
			}
			else
			{
				pMachines[i] = 3;
			}
		}
	}
	if (pType == 2)
	{
		for (int i = 0; i < pStage; i++)
		{
			if (i == 1)
			{
				pMachines[i] = 1;
			}
			else
			{
				pMachines[i] = 3;
			}
		}
	}
	if (pType == 3)
	{
		for (int i = 0; i < pStage; i++)
		{
			if (i == 1)
			{
				pMachines[i] = 2;
			}
			else
			{
				pMachines[i] = 3;
			}
		}
	}
	if (pType == 4)
	{
		for (int i = 0; i < pStage; i++)
		{
			pMachines[i] = 3;
		}
	}*/
	machineSpeed.resize(InStage);
	for (int i = 0; i < InStage; i++)
	{
		machineSpeed[i].resize(5);    //number of speed levels: 5  machineSpeed[i] = 1 + rand() % 5
	}

	for (int i = 0; i < InStage; i++)
	{
		bool isOK = false;
		while (!isOK)
		{
			for (int v = 0; v < 5; v++)
			{
				//machineSpeed[i][v] = (int)((rand() * 1.0 / RAND_MAX + 1.0) * 100) * 1.0 / 100;
				machineSpeed[i][v] = 1 + rand() % 5;
				//cout << machineSpeed[i][v] << endl;

			}
			for (int v = 1; v < 5; v++)
			{
				if (machineSpeed[i][v] <= machineSpeed[i][v - 1])
				{
					isOK = false;
					break;
				}
				isOK = true;
			}
		}
	}

	/*pUnitTime.resize(InStage);
	for (int i = 0; i < pStage; i++)
	{
		pUnitTime[i].resize(pJob);
	}
	for (int i = 0; i < pStage; i++)
	{
		for (int j = 0; j < pJob; j++)
		{

			pUnitTime[i][j] = 1 + rand() % 10;
		}
	}*/
	bTime.resize(InStage);
	for (int i = 0; i < InStage; i++)
	{
		bTime[i].resize(InJob);
	}
	for (int i = 0; i < InStage; i++)
	{
		for (int j = 0; j < InJob; j++)
		{
			bTime[i][j] = 1 + rand() % 10;  //[1-99]
		}
	}
	
	pUnitTime.resize(pStage);
	for (int i = 0; i < pStage; i++)
	{
		pUnitTime[i].resize(pJob);
	}
	for (int i = 0; i < InStage; i++)
	{
		for (int j = 0; j < InJob; j++)
		{
			pUnitTime[i][j].resize(5);    //five speed levels
		}
	}
	for (int k = 0; k < pStage; k++)
	{
			for (int j = 0; j < pJob; j++)

			{
				for (int v = 0; v < 5; v++)
				{
					pUnitTime[k][j][v] = (int)(bTime[k][j]/ machineSpeed[k][v] * 100) * 1.0 / 100;
					//cout << pUnitTime[k][j][v] << endl;
				}
			}
	}

	pLot.resize(pJob);
	for (int j = 0; j < pJob; j++)
	{
		pLot[j] = 50 + rand() % 51;
	}
	pSetupTime.resize(pStage);
	for (int i = 0; i < pStage; i++)
	{
		pSetupTime[i].resize(pJob);
	}
	for (int i = 0; i < pStage; i++)
	{
		for (int j = 0; j < pJob; j++)
		{

			pSetupTime[i][j] = 50 + rand() % 51;

		}
	}
	pTransferTime.resize(pStage);
	for (int i = 0; i < pStage; i++)
	{
		pTransferTime[i].resize(pJob);
	}

	for (int i = 0; i < pStage; i++)
	{
		for (int j = 0; j < pJob; j++)
		{
			pTransferTime[i][j] = 5 + rand() % 11;
		}
	}
	pPower.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		pPower[k].resize(pJob);
	}
	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{
			pPower[k][j].resize(5);
		}
	}
	for (int k = 0; k < pStage; k++)
	{
			for (int j = 0; j < pJob; j++)
			{
				for (int v = 0; v < 5; v++)
				{
					pPower[k][j][v] = (int)(4 * machineSpeed[k][v] * machineSpeed[k][v] * 100) * 1.0 / 100;
				}
			}
	}
	ip = 1;
	sp = 2;
}



void OutputInstances(int InJob, int InStage, int InType, int NumberOfInstance)
{
	char InsFile[80];
	sprintf_s(InsFile, "Instances\\data%d%d%d.txt", NumberOfInstance);

	fstream fout(InsFile, fstream::out | fstream::app);



	fout << MaxMS << " " << MinMS << " " << MaxTEC << " " << MinTEC << endl;

	fout << pJob << " " << pStage << endl;


	for (int i = 0; i < pStage; i++)
	{
		fout << pMachines[i] << " ";
	}
	fout << endl;



	for (int i = 0; i < pStage; i++)
	{
		for (int j = 0; j < pJob; j++)
		{
			for (int v = 0; v < 5; v++)
			{
				fout << pUnitTime[i][j][v] << " ";
			}
		}
		fout << endl;
	}

	for (int i = 0; i < pStage; i++)
	{
		for (int j = 0; j < pJob; j++)
		{
			fout << bTime[i][j] << "";  //[1-99]
		}
		fout << endl;
	}

	for (int j = 0; j < pJob; j++)
	{
		fout << pLot[j] << " ";
	}
	fout << endl;


	for (int i = 0; i < pStage; i++)
	{
		for (int j = 0; j < pJob; j++)
		{

			fout << pSetupTime[i][j] << " ";

		}
		fout << endl;
	}


	for (int i = 0; i < pStage; i++)
	{
		for (int j = 0; j < pJob; j++)
		{

			fout << pTransferTime[i][j] << " ";

		}
		fout << endl;
	}


	fout.close();
}

void GetUpperAndLowerBounds(int InType)
{
	//It is natrual to determien the minimum and maximum totoal number of sublots
	MinTEC = 0.0;
	MaxTEC = 0.0;
	MinMS = 0;
	MaxMS = 0;
	//get the MinMS  with extreme lot streaming 
	vector<vector<vector<int>>> IndependentMS;
	IndependentMS.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		IndependentMS[k].resize(pJob);
	}
	for (int k = 0; k < pStage; k++)
	{
		for (int i = 0; i < pJob; i++)
		{
			IndependentMS[k][i].resize(pLot[i], 0);
		}
	}

	//for the first stage
	for (int i = 0; i < pJob; i++)
	{
		for (int j = 0; j < pLot[i]; j++)
		{

			for (int v = 0; v < 5; v++)
			{
				if (j == 0)
				{
					IndependentMS[0][i][j] = pSetupTime[0][i] + pUnitTime[0][i][4];
				}
				else
				{
					IndependentMS[0][i][j] = IndependentMS[0][i][j - 1] + pUnitTime[0][i][4];
				}

			}

		}
	}

	//for the other stages
	for (int k = 1; k < pStage; k++)
	{
		for (int i = 0; i < pJob; i++)
		{
			for (int j = 0; j < pLot[i]; j++)
			{
				for (int v = 0; v < 5; v++)
				{
					if (j == 0)
					{
						if (pSetupTime[k][i] > IndependentMS[k - 1][i][j] + pTransferTime[k][i])
						{
							IndependentMS[k][i][j] = pSetupTime[k][i] + pUnitTime[k][i][4];
						}
						else
						{
							IndependentMS[k][i][j] = IndependentMS[k - 1][i][j] + pTransferTime[k][i] + pUnitTime[k][i][4];
						}
					}
					else
					{
						IndependentMS[k][i][j] = max(IndependentMS[k][i][j - 1], IndependentMS[k - 1][i][j] + pTransferTime[k][i]) + pUnitTime[k][i][4];
					}
				}
			}
		}
	}


	MinMS = IndependentMS[pStage - 1][0][pLot[0] - 1];

	for (int i = 1; i < pJob; i++)
	{
		if (IndependentMS[pStage - 1][i][pLot[i] - 1] > MinMS)
		{
			MinMS = IndependentMS[pStage - 1][i][pLot[i] - 1];
		}
	}


	//get the MaxMS with no lot streaming
	vector<vector<int>> IndependentMS1;
	IndependentMS1.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		IndependentMS1[k].resize(pJob, 0);
	}


	//for the first stage
	for (int i = 0; i < pJob; i++)
	{
		for (int e = 0; e < MaxSublotQuantity; e++)
		{
			for (int k = 0; k < 5; k++)
			{

				IndependentMS1[0][i] = pSetupTime[0][i] + pUnitTime[0][i][0] * pLot[i];
			}
		}
	}

	//for the other stages
	for (int k = 1; k < pStage; k++)
	{
		for (int i = 0; i < pJob; i++)
		{

			if (pSetupTime[k][i] > IndependentMS1[k - 1][i] + pTransferTime[k][i])
			{
				IndependentMS1[k][i] = pSetupTime[k][i] + pUnitTime[k][i][0] * pLot[i];
			}
			else
			{
				IndependentMS1[k][i] = IndependentMS1[k - 1][i] + pUnitTime[k][i][0] * pLot[i] + pTransferTime[k][i];
			}
		}
	}


	for (int i = 0; i < pJob; i++)
	{
		MaxMS += IndependentMS1[pStage - 1][i];
	}
	//get lowtec
	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{

			MinTEC += pUnitTime[k][j][0] * pLot[j] * pPower[k][j][0];

		}
	}
	//get MAXTEC
	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{
			MaxTEC += pUnitTime[k][j][4] * pLot[j] * pPower[k][j][4];
			MaxTEC += 2 * pSetupTime[k][j];
		}
	}
}
class Pair
{
public:
	int dim;
	int value[MaxSublotQuantity];
};


//子批同时到达
class PairLess
{
public:
	bool operator() (Pair& a, Pair& b)
	{
		bool result = false;
		int i = 0;
		do
		{

			if (a.value[i] < b.value[i])
			{
				result = true;
			}
			else if (a.value[i] > b.value[i])
			{
				result = false;
			}
			i++;
		} while (a.value[i - 1] == b.value[i - 1] && i < MaxSublotQuantity);

		return result;

	}
};
//批次同时到达
class PairLess2
{
public:
	bool operator() (Pair& a, Pair& b)
	{
		if (a.value[0] < b.value[0])
		{
			return true;
		}
		else if (a.value[0] == b.value[0])
		{
			if (a.value[1] < b.value[1])
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}

	}
};


//不同时到达
class PairLess3
{
	bool operator()(Pair& a, Pair& b)
	{
		return a.value[0] < b.value[0];
	}
};


class PairX
{
public:
	int dim;
	int value;
};

class PairXLess
{
public:
	bool operator () (PairX& a1, PairX& a2)
	{
		return a1.value < a2.value;
	}
};

class PairXGreater
{
public:
	bool operator () (PairX& a1, PairX& a2)
	{
		return a1.value > a2.value;
	}
};
