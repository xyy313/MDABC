#pragma once
#include <vector>
#include "HFSP_CV.h"
#include <algorithm>
#include "Global.h"

using namespace std;

class Individual
{
public:
	//the solution representation
	vector<int> pJobSeq;  //job
	vector<vector<int>> pJobSplit;  //job-sublot
	vector<vector<vector<int>>> pSpeedSelection;//stage-job-sublot
	int pAge;
	int neighborhoodAge;
public:
	int MS;				   //the objective: makespan
	int TEC;				   //the objective: total energy cost objective
	double Fitness;            //the fitness for the MOEA/D
	int n;					   //记录支配this的个体数目
	double distance;			//crowding distance
	double topsis_pdistance;      //topsis_distance
	double topsis_ndistance;      //topsis_ndistance
	double topsis_rdistance;       //relative distance

	double WeightMS;
	double WeightTEC;

public:
	vector<int> IndexNeighbor; //neighbor index
	int CurrentWeightIndex;//当前的权重索引
	vector<double> NeighborWeightMS;  //neighbor weight MS
	vector<double> NeighborWeightTEC; //neighbhor weight TEC

	int CurrentNeighbrohoodIndex;     //current neighborhood structure index

public:
	int NS;    // the index of neighborhood structure
	int NSCF;   //consectve failure to update the current solution

public:
	Individual();    //constructor
	Individual(int k);
	~Individual();    //deconstructor

	void getObjectives2(int DecodingMethod);       //obtain the objectives  机器选择规则：最先完成 （这个要比1好）
	void getObjectives4(int DecodingMethod);       //obtain the objectives 对objective3()进行修正 修正成功  综上所述 2和4是正确的
	void getFitnessValue();
public:
	void conductPermution();

public:
	bool operator>(Individual& obj);   //dominate
	bool operator==(Individual& obj);  //equal
	bool operator<(Individual& obj);   //dominated
	//Individual & operator=(Individual &obj);
public:
	void getNeigbors(Individual& Neighbor, int k);    //generate a neighbor
};

Individual::Individual()
{

}

Individual::Individual(int k)
{

	if (k == 1)
	{
		//Initialzie the pJobSeq
		pJobSeq.resize(pJob);
		for (int i = 0; i < pJob; i++)
		{
			pJobSeq[i] = i;
		}
		random_shuffle(pJobSeq.begin(), pJobSeq.end());

		//Initialize the pJobSplit
		//产生 pSublotJobs
		pJobSplit.resize(pJob);
		for (int i = 0; i < pJob; i++)
		{
			pJobSplit[i].resize(MaxSublotQuantity);
		}
		for (int i = 0; i < pJob; i++)
		{
			for (int j = 0; j < MaxSublotQuantity; j++)
			{
				pJobSplit[i][j] = floor(pLot[i] / MaxSublotQuantity);
			}
		}
		vector<int> r;
		r.resize(pJob);
		for (int i = 0; i < pJob; i++)
		{
			int sum = 0;
			for (int j = 0; j < MaxSublotQuantity; j++)
			{
				sum += pJobSplit[i][j];
			}
			r[i] = pLot[i] - sum;

			//随机选择一个子批加上r[i]
			int rIndex = rand() % MaxSublotQuantity;
			pJobSplit[i][rIndex] += r[i];
		}
		////initialize the speed selection matrix
		pSpeedSelection.resize(pStage);
		for (int k = 0; k < pStage; k++)
		{
			pSpeedSelection[k].resize(pJob);
		}
		for (int k = 0; k < pStage; k++)
		{
			for (int j = 0; j < pJob; j++)
			{
				pSpeedSelection[k][j].resize(MaxSublotQuantity);
			}
		}

		for (int k = 0; k < pStage; k++)
		{
			for (int j = 0; j < pJob; j++)
			{
				for (int e = 0; e < MaxSublotQuantity; e++)
				{
					pSpeedSelection[k][j][e] = rand() % 5;
				}
			}
		}
		pAge = 0;
		neighborhoodAge = 0;
		CurrentWeightIndex = 0;
		CurrentNeighbrohoodIndex = 1;
		getObjectives2(1);
	}
	if (k == 2)
	{
		//Initialzie the pJobSeq
		pJobSeq.resize(pJob);
		for (int i = 0; i < pJob; i++)
		{
			pJobSeq[i] = i;
		}
		random_shuffle(pJobSeq.begin(), pJobSeq.end());

		//Initialize the pJobSplit
		//产生 pSublotJobs
		pJobSplit.resize(pJob);
		for (int i = 0; i < pJob; i++)
		{
			pJobSplit[i].resize(MaxSublotQuantity);
		}
		for (int i = 0; i < pJob; i++)
		{
			for (int j = 0; j < MaxSublotQuantity; j++)
			{
				pJobSplit[i][j] = floor(pLot[i] / MaxSublotQuantity);
			}
		}
		vector<int> r;
		r.resize(pJob);
		for (int i = 0; i < pJob; i++)
		{
			int sum = 0;
			for (int j = 0; j < MaxSublotQuantity; j++)
			{
				sum += pJobSplit[i][j];
			}
			r[i] = pLot[i] - sum;

			//随机选择一个子批加上r[i]
			int rIndex = rand() % MaxSublotQuantity;
			pJobSplit[i][rIndex] += r[i];
		}
		////initialize the speed selection matrix
		pSpeedSelection.resize(pStage);
		for (int k = 0; k < pStage; k++)
		{
			pSpeedSelection[k].resize(pJob);
		}
		for (int k = 0; k < pStage; k++)
		{
			for (int j = 0; j < pJob; j++)
			{
				pSpeedSelection[k][j].resize(MaxSublotQuantity);
			}
		}
		for (int k = 0; k < pStage; k++)
		{
			for (int j = 0; j < pJob; j++)
			{
				for (int e = 0; e < MaxSublotQuantity; e++)
				{
					pSpeedSelection[k][j][e] = rand() % 5;
				}
			}
		}
		pAge = 0;
		neighborhoodAge = 0;
		CurrentWeightIndex = 0;
		CurrentNeighbrohoodIndex = 1;
		getObjectives4(2);
	}
}
//l邻域结构的设计
void Individual::getNeigbors(Individual& Neighbor, int k)
{
	bool IsImprovement = false;
	//int number = rand() % 3 + 1;   //3的值需要做实验就行设定
	int number = 1;
	Individual tempIndividual = *this;
	//Individual *tempIndividual = new Individual(*this);

	tempIndividual.pAge = 0;

	//job sequence insert
	if (k == 3)
	{
		for (int i = 0; i < number; i++)
		{
			int a = rand() % pJob;
			int b = rand() % pJob;
			while (a == b)
			{
				a = rand() % pJob;
				b = rand() % pJob;
			}

			if (b > a)
			{
				int tempValue = tempIndividual.pJobSeq[b];
				for (int i = b; i > a; i--)
				{
					tempIndividual.pJobSeq[i] = tempIndividual.pJobSeq[i - 1];
				}
				tempIndividual.pJobSeq[a] = tempValue;
			}
			else
			{
				int tempValue = tempIndividual.pJobSeq[b];
				for (int i = b; i < a; i++)
				{
					tempIndividual.pJobSeq[i] = tempIndividual.pJobSeq[i + 1];
				}
				tempIndividual.pJobSeq[a] = tempValue;
			}
		}
	}
	//job sequence swap
	else if (k == 4)
	{
		for (int i = 0; i < number; i++)
		{
			int a = rand() % pJob;
			int b = rand() % pJob;
			while (a == b)
			{
				a = rand() % pJob;
				b = rand() % pJob;
			}

			swap(tempIndividual.pJobSeq[a], tempIndividual.pJobSeq[b]);
		}
	}
	//批次分割
	else if (k == 5)
	{
		for (int i = 0; i < number; i++)
		{
			//任选一个批次
			int r = rand() % pJob;

			//任选这个批次的两个子批
			int a, b;
			do
			{
				a = rand() % MaxSublotQuantity;
				b = rand() % MaxSublotQuantity;
			} while (a == b);

			if (tempIndividual.pJobSplit[r][a] > 0 || tempIndividual.pJobSplit[r][b] > 0)
			{
				//定义一个步长
				int step = rand() % 3 + 1;
				double r2 = rand() * 1.0 / RAND_MAX;
				if (r2 < 0.5)
				{
					tempIndividual.pJobSplit[r][a] += step;
					tempIndividual.pJobSplit[r][b] -= step;
				}
				else
				{
					tempIndividual.pJobSplit[r][a] += step;
					tempIndividual.pJobSplit[r][b] -= step;
				}

				if (tempIndividual.pJobSplit[r][0] == 0 && tempIndividual.pJobSplit[r][a] < 0 && tempIndividual.pJobSplit[r][b] < 0)
				{
					tempIndividual.pJobSplit[r][a] = this->pJobSplit[r][a];
					tempIndividual.pJobSplit[r][b] = this->pJobSplit[r][b];
				}
			}
		}
	}
	//job insert+job split
	else if (k == 1)
	{
		for (int i = 0; i < number; i++)
		{
			int a = rand() % pJob;
			int b = rand() % pJob;
			while (a == b)
			{
				a = rand() % pJob;
				b = rand() % pJob;
			}
			if (b > a)
			{
				int tempValue = tempIndividual.pJobSeq[b];
				for (int i = b; i > a; i--)
				{
					tempIndividual.pJobSeq[i] = tempIndividual.pJobSeq[i - 1];
				}
				tempIndividual.pJobSeq[a] = tempValue;
			}
			else
			{
				int tempValue = tempIndividual.pJobSeq[b];
				for (int i = b; i < a; i++)
				{
					tempIndividual.pJobSeq[i] = tempIndividual.pJobSeq[i + 1];
				}
				tempIndividual.pJobSeq[a] = tempValue;
			}
			for (int k = 0; k < pStage; k++)
			{
				for (int j = 0; j < pJob; j++)
				{
					int tempstage = rand() % pStage;
					int tempjob = rand() % pJob;
				}
			}
			//任选一个子批数量大于1的批次
			int r = rand() % pJob;

			//任选这个批次的两个子批
			do
			{
				a = rand() % MaxSublotQuantity;
				b = rand() % MaxSublotQuantity;
			} while (a == b);
			if (tempIndividual.pJobSplit[r][a] > 0 || tempIndividual.pJobSplit[r][b] > 0)
			{
				//定义一个步长
				int step = rand() % 3 + 1;
				double r2 = rand() * 1.0 / RAND_MAX;
				if (r2 < 0.5)
				{
					tempIndividual.pJobSplit[r][a] += step;
					tempIndividual.pJobSplit[r][b] -= step;
				}
				else
				{
					tempIndividual.pJobSplit[r][a] -= step;
					tempIndividual.pJobSplit[r][b] += step;
				}

				if (tempIndividual.pJobSplit[r][0] == 0 && tempIndividual.pJobSplit[r][a] < 0 && tempIndividual.pJobSplit[r][b] < 0)
				{
					tempIndividual.pJobSplit[r][a] = this->pJobSplit[r][a];
					tempIndividual.pJobSplit[r][b] = this->pJobSplit[r][b];
				}
			}
		}
		tempIndividual.getObjectives2(1);
		tempIndividual.getFitnessValue();
	}
	//job swap +job split
	else if (k == 2)
	{
		for (int i = 0; i < number; i++)
		{
			int a = rand() % pJob;
			int b = rand() % pJob;
			while (a == b)
			{
				a = rand() % pJob;
				b = rand() % pJob;
			}

			swap(tempIndividual.pJobSeq[a], tempIndividual.pJobSeq[b]);
			for (int k = 0; k < pStage; k++)
			{
				for (int j = 0; j < pJob; j++)
				{
					int tempstage = rand() % pStage;
					int tempjob = rand() % pJob;

				}
			}
			//任选一个子批数量大于1的批次
			int r = rand() % pJob;

			//任选这个批次的两个子批
			do
			{
				a = rand() % MaxSublotQuantity;
				b = rand() % MaxSublotQuantity;
			} while (a == b);
			if (tempIndividual.pJobSplit[r][a] > 0 || tempIndividual.pJobSplit[r][b] > 0)
			{
				//定义一个步长
				int step = rand() % 3 + 1;
				double r2 = rand() * 1.0 / RAND_MAX;
				if (r2 < 0.5)
				{
					tempIndividual.pJobSplit[r][a] += step;
					tempIndividual.pJobSplit[r][b] -= step;
				}
				else
				{
					tempIndividual.pJobSplit[r][a] -= step;
					tempIndividual.pJobSplit[r][b] += step;
				}

				if (tempIndividual.pJobSplit[r][0] == 0 && tempIndividual.pJobSplit[r][a] < 0 && tempIndividual.pJobSplit[r][b] < 0)
				{
					tempIndividual.pJobSplit[r][a] = this->pJobSplit[r][a];
					tempIndividual.pJobSplit[r][b] = this->pJobSplit[r][b];
				}

			}
		}
		tempIndividual.getObjectives2(1);
		tempIndividual.getFitnessValue();
	}
	//speedselection
	else  if (k == 6)
	{
		for (int i = 0; i < number; i++)
		{
			int tempstage = rand() % pStage;
			int tempjob = rand() % pJob;
			int tempe = rand() % MaxSublotQuantity;
			tempIndividual.pSpeedSelection[tempstage][tempjob][tempe] = rand() % 5;
		}
		tempIndividual.getObjectives2(1);
		tempIndividual.getFitnessValue();
	}
	//job swap +speed selection
	else  if (k == 7)
	{
		for (int i = 0; i < number; i++)
		{
			int a = rand() % pJob;
			int b = rand() % pJob;
			while (a == b)
			{
				a = rand() % pJob;
				b = rand() % pJob;
			}
			swap(tempIndividual.pJobSeq[a], tempIndividual.pJobSeq[b]);
			for (int i = 0; i < number; i++)
			{
				int tempstage = rand() % pStage;
				int tempjob = rand() % pJob;
				int tempe = rand() % MaxSublotQuantity;
				tempIndividual.pSpeedSelection[tempstage][tempjob][tempe] = rand() % 5;

			}
		}
		tempIndividual.getObjectives2(1);
		tempIndividual.getFitnessValue();
	}
	else  if (k == 8)
	{
		for (int i = 0; i < number; i++)
		{
			int a = rand() % pJob;
			int b = rand() % pJob;
			while (a == b)
			{
				a = rand() % pJob;
				b = rand() % pJob;
			}

			if (b > a)
			{
				int tempValue = tempIndividual.pJobSeq[b];
				for (int i = b; i > a; i--)
				{
					tempIndividual.pJobSeq[i] = tempIndividual.pJobSeq[i - 1];
				}
				tempIndividual.pJobSeq[a] = tempValue;
			}
			else
			{
				int tempValue = tempIndividual.pJobSeq[b];
				for (int i = b; i < a; i++)
				{
					tempIndividual.pJobSeq[i] = tempIndividual.pJobSeq[i + 1];
				}
				tempIndividual.pJobSeq[a] = tempValue;
			}
			for (int i = 0; i < number; i++)
			{
				int tempstage = rand() % pStage;
				int tempjob = rand() % pJob;
				int tempe = rand() % MaxSublotQuantity;
				tempIndividual.pSpeedSelection[tempstage][tempjob][tempe] = rand() % 5;
			}
		}
		tempIndividual.getObjectives2(1);
		tempIndividual.getFitnessValue();
	}
}
void Individual::getObjectives2(int DecodingMethod)  //obtain the objectives
{
	vector<int> Seq = pJobSeq;

	vector<vector<vector<int>>> machineAssignment;  //机器分配
	machineAssignment.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		machineAssignment[k].resize(pMachines[k]);
	}

	vector<vector<int>> mIdleTime;//找出最早空闲的机器
	mIdleTime.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		mIdleTime[k].resize(pMachines[k]);
	}
	//all machines are availale at 0  其实可以不用要 默认都为0
	for (int k = 0; k < pStage; k++)
	{
		for (int i = 0; i < pMachines[k]; i++)
		{
			mIdleTime[k][i] = 0;
		}
	}
	vector<vector<vector<int>>> STime, CTime;
	STime.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		STime[k].resize(pJob);
	}

	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{
			STime[k][j].resize(MaxSublotQuantity);
		}
	}
	CTime = STime;
	//Schedule the first stage
	for (int j = 0; j < pJob; j++)
	{
		//Select a machine
		int mpt;
		int	minIdle = INT_MAX;
		for (int i = 0; i < pMachines[0]; i++)
		{
			if (mIdleTime[0][i] < minIdle)
			{
				minIdle = mIdleTime[0][i];
				mpt = i;
			}
		}
		machineAssignment[0][mpt].push_back(Seq[j]);//
		//工件批量的调度
		for (int e = 0; e < MaxSublotQuantity; e++)
		{
			//第一个子批有启动时间
			if (e == 0)
			{
				STime[0][Seq[j]][e] = mIdleTime[0][mpt] + pSetupTime[0][Seq[j]];//序列不相关-》序列相关
				CTime[0][Seq[j]][e] = STime[0][Seq[j]][e] + pUnitTime[0][Seq[j]][pSpeedSelection[0][Seq[j]][e]] * pJobSplit[Seq[j]][e];
				mIdleTime[0][mpt] = CTime[0][Seq[j]][e];

			}
			else
			{
				STime[0][Seq[j]][e] = mIdleTime[0][mpt];
				CTime[0][Seq[j]][e] = STime[0][Seq[j]][e] + pUnitTime[0][Seq[j]][pSpeedSelection[0][Seq[j]][e]] * pJobSplit[Seq[j]][e];
				mIdleTime[0][mpt] = CTime[0][Seq[j]][e];

			}

		}

	}

	//Schedule jobs in other stages as earlier as possible
	for (int k = 1; k < pStage; k++)
	{
		vector<Pair> ch;
		ch.resize(pJob);
		//Pair<double>* ch = new Pair<double>[pJob];


		if (DecodingType == 1)
		{
			//Sorting jobs according to the completion time of the first sublot of jobs in previous stage
			for (int j = 0; j < pJob; j++)
			{
				ch[j].dim = j;
				for (int e = 0; e < MaxSublotQuantity; e++)
				{
					ch[j].value[e] = CTime[k - 1][j][e];  //子批优先
				}

			}

			sort(ch.begin(), ch.end(), PairLess());
			for (int j = 0; j < pJob; j++)
			{
				Seq[j] = ch[j].dim;
			}
		}
		else if (DecodingType == 2)
		{
			//Sorting jobs according to the completion time of the first sublot of jobs in previous stage
			for (int j = 0; j < pJob; j++)
			{
				for (int v = 0; v < 5; v++)
				{
					ch[j].dim = j;
					ch[j].value[0] = CTime[k - 1][j][MaxSublotQuantity - 1];  //批次j的最后一个子批在上一阶段的完成时间,批次优先
					if (pUnitTime[k - 1][j][v] > pUnitTime[k][j][v])
					{
						ch[j].value[1] = pUnitTime[k - 1][j][v] - pUnitTime[k][j][v];
					}
					else
					{
						ch[j].value[1] = 0;
					}
				}
			}
			sort(ch.begin(), ch.end(), PairLess2());
			for (int j = 0; j < pJob; j++)
			{
				Seq[j] = ch[j].dim;
			}
		}

		//schedule jobs
		for (int j = 0; j < pJob; j++)
		{
			//Select a machine
			int mpt;

			int minIdle = INT_MAX;
			for (int i = 0; i < pMachines[k]; i++)
			{
				if (mIdleTime[k][i] < minIdle)
				{
					minIdle = mIdleTime[k][i];
					mpt = i;
				}
			}

			machineAssignment[k][mpt].push_back(Seq[j]);
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				//CTime[k][JobSeq[j]] = max(mIdleTime[k][machineIndex] + setupTime[k][JobSeq[j]][JobSeq[j]], CTime[k - 1][JobSeq[j]]) + bTime[k][JobSeq[j]];

				if (e == 0)
				{
					CTime[k][Seq[j]][e] = max(mIdleTime[k][mpt] + pSetupTime[k][Seq[j]], CTime[k - 1][Seq[j]][e] + pTransferTime[k - 1][Seq[j]]) + pUnitTime[k][Seq[j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
					STime[k][Seq[j]][e] = CTime[k][Seq[j]][e] - pUnitTime[k][Seq[j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
					mIdleTime[k][mpt] = CTime[k][Seq[j]][e];
				}
				else
				{
					STime[k][Seq[j]][e] = max(mIdleTime[k][mpt], CTime[k - 1][Seq[j]][e] + pTransferTime[k - 1][Seq[j]] * ((pJobSplit[Seq[j]][e] > 0) ? 1 : 0));
					CTime[k][Seq[j]][e] = STime[k][Seq[j]][e] + pUnitTime[k][Seq[j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
					mIdleTime[k][mpt] = CTime[k][Seq[j]][e];
				}

			}

			//get makespan
			MS = 0;

			for (int i = 0; i < pJob; i++)
			{
				if (CTime[pStage - 1][i][MaxSublotQuantity - 1] > MS)
				{
					MS = CTime[pStage - 1][i][MaxSublotQuantity - 1];
				}

			}
			//get energyconsumption
			TEC = 0.0;
			double EC1 = 0.0;   //energy cost in the processing mode
			double EC2 = 0.0;   //energy cost in the setup mode
			double EC3 = 0.0;   //energy cost in the standby mode

			for (int k = 0; k < pStage; k++)
			{
				for (int i = 0; i < pMachines[k]; i++)
				{
					double totalPT = 0.0;
					for (int j = 0; j < machineAssignment[k][i].size(); j++)
					{
						for (int e = 0; e < MaxSublotQuantity; e++)
						{
							EC1 += pUnitTime[k][machineAssignment[k][i][j]][pSpeedSelection[k][Seq[j]][e]] * pPower[k][machineAssignment[k][i][j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
							totalPT += pUnitTime[k][machineAssignment[k][i][j]][pSpeedSelection[k][Seq[j]][e]];
							EC2 += pSetupTime[k][machineAssignment[k][i][j]] * sp;
							totalPT += pSetupTime[k][machineAssignment[k][i][j]];
						}
					}
					double IdleTime = mIdleTime[k][i] - totalPT;
					EC3 += IdleTime * ip;
				}
			}
			TEC = EC1 + EC2 + EC3;
		}
	}
}
void Individual::getObjectives4(int DecodingMethod)  //obtain the objectives
{
	vector<int> Seq = pJobSeq;
	vector<vector<vector<int>>> machineAssignment;  //机器分配
	machineAssignment.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		machineAssignment[k].resize(pMachines[k]);
	}
	vector<vector<int>> mIdleTime;//找出最早空闲的机器
	mIdleTime.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		mIdleTime[k].resize(pMachines[k]);
	}
	//all machines are availale at 0  其实可以不用要 默认都为0
	for (int k = 0; k < pStage; k++)
	{
		for (int i = 0; i < pMachines[k]; i++)
		{
			mIdleTime[k][i] = 0;
		}
	}
	vector<vector<vector<int>>> STime, CTime, TSTime;
	STime.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		STime[k].resize(pJob);
	}

	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{
			STime[k][j].resize(MaxSublotQuantity);
		}
	}
	CTime = STime;
	TSTime = STime;
	//Schedule the first stage
	for (int j = 0; j < pJob; j++)
	{
		//Select a machine
		int mpt;
		int	minIdle = INT_MAX;
		for (int i = 0; i < pMachines[0]; i++)
		{
			if (mIdleTime[0][i] < minIdle)
			{
				minIdle = mIdleTime[0][i];
				mpt = i;
			}
		}
		machineAssignment[0][mpt].push_back(Seq[j]);//
		//工件批量的调度
		for (int e = 0; e < MaxSublotQuantity; e++)
		{
			//第一个子批有启动时间
			if (e == 0)
			{
				STime[0][Seq[j]][e] = mIdleTime[0][mpt] + pSetupTime[0][Seq[j]];
				CTime[0][Seq[j]][e] = STime[0][Seq[j]][e] + pUnitTime[0][Seq[j]][pSpeedSelection[0][Seq[j]][e]] * pJobSplit[Seq[j]][e];
				mIdleTime[0][mpt] = CTime[0][Seq[j]][e];


			}
			else
			{
				STime[0][Seq[j]][e] = mIdleTime[0][mpt];
				CTime[0][Seq[j]][e] = STime[0][Seq[j]][e] + pUnitTime[0][Seq[j]][pSpeedSelection[0][Seq[j]][e]] * pJobSplit[Seq[j]][e];
				mIdleTime[0][mpt] = CTime[0][Seq[j]][e];

			}
			TSTime[0][Seq[j]][e] = mIdleTime[0][mpt] - pUnitTime[0][Seq[j]][pSpeedSelection[0][Seq[j]][e]] * pJobSplit[Seq[j]][e];
		}

	}
	//for the other stages
	for (int k = 1; k < pStage; k++)
	{
		vector<Pair> ch;
		ch.resize(pJob);
		if (DecodingType == 1)
		{
			//Sorting jobs according to the completion time of the first sublot of jobs in previous stage
			for (int j = 0; j < pJob; j++)
			{
				ch[j].dim = j;
				for (int e = 0; e < MaxSublotQuantity; e++)
				{
					ch[j].value[e] = CTime[k - 1][j][e];  //子批优先
				}

			}
			sort(ch.begin(), ch.end(), PairLess());
			for (int j = 0; j < pJob; j++)
			{
				Seq[j] = ch[j].dim;
			}
		}
		else if (DecodingType == 2)
		{
			//Sorting jobs according to the completion time of the first sublot of jobs in previous stage
			for (int j = 0; j < pJob; j++)
			{
				for (int v = 0; v < 5; v++)
				{
					ch[j].dim = j;
					ch[j].value[0] = CTime[k - 1][j][MaxSublotQuantity - 1];  //批次j的最后一个子批在上一阶段的完成时间,批次优先
					if (pUnitTime[k - 1][j][v] > pUnitTime[k][j][v])
					{
						ch[j].value[1] = pUnitTime[k - 1][j][v] - pUnitTime[k][j][v];
					}
					else
					{
						ch[j].value[1] = 0;
					}
				}
			}
			sort(ch.begin(), ch.end(), PairLess2());
			for (int j = 0; j < pJob; j++)
			{
				Seq[j] = ch[j].dim;
			}
		}

		//schedule jobs
		for (int j = 0; j < pJob; j++)
		{
			//Select a machine
			int mpt;

			int minIdle = INT_MAX;
			for (int i = 0; i < pMachines[k]; i++)
			{
				if (mIdleTime[k][i] < minIdle)
				{
					minIdle = mIdleTime[k][i];
					mpt = i;
				}
			}

			machineAssignment[k][mpt].push_back(Seq[j]);
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				//CTime[k][JobSeq[j]] = max(mIdleTime[k][machineIndex] + setupTime[k][JobSeq[j]][JobSeq[j]], CTime[k - 1][JobSeq[j]]) + bTime[k][JobSeq[j]];
				if (e == 0)
				{
					CTime[k][Seq[j]][e] = max(mIdleTime[k][mpt] + pSetupTime[k][Seq[j]], CTime[k - 1][Seq[j]][e] + pTransferTime[k - 1][Seq[j]]) + pUnitTime[k][Seq[j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
					STime[k][Seq[j]][e] = CTime[k][Seq[j]][e] - pUnitTime[k][Seq[j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
					mIdleTime[k][mpt] = CTime[k][Seq[j]][e];
				}
				else
				{
					STime[k][Seq[j]][e] = max(mIdleTime[k][mpt], CTime[k - 1][Seq[j]][e] + pTransferTime[k - 1][Seq[j]] * ((pJobSplit[Seq[j]][e] > 0) ? 1 : 0));
					CTime[k][Seq[j]][e] = STime[k][Seq[j]][e] + pUnitTime[k][Seq[j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
					mIdleTime[k][mpt] = CTime[k][Seq[j]][e];
				}
				TSTime[k][Seq[j]][e] = STime[k][Seq[j]][e] + pSetupTime[k][Seq[j]];
			}
			//get makespan
			MS = 0;

			for (int i = 0; i < pJob; i++)
			{
				if (CTime[pStage - 1][i][MaxSublotQuantity - 1] > MS)
				{
					MS = CTime[pStage - 1][i][MaxSublotQuantity - 1];
				}

			}
			//backward decoding 
		   //for the last stage
			for (int i = 0; i < pMachines[pStage - 1]; i++)
			{
				for (int j = machineAssignment[pStage - 1][i].size() - 1; j >= 0; j--)
				{
					for (int e = 0; e < MaxSublotQuantity; e++)
					{
						int NextJobSTime;
						if (j == machineAssignment[pStage - 1][i].size() - 1)
						{
							NextJobSTime = INT_MAX;
						}
						else
						{
							NextJobSTime = STime[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][e] * pJobSplit[Seq[j]][e];
						}
						int speedLevel = pSpeedSelection[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][e];
						for (int v = 0; v < speedLevel; v++)
						{
							double difference;
							difference = pUnitTime[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][v] * pJobSplit[Seq[j]][e] - pUnitTime[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][speedLevel] * pJobSplit[Seq[j]][e];
							double NovelCTime = CTime[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][e] * pJobSplit[Seq[j]][e] + difference;
							if (NovelCTime <= min(NextJobSTime, MS))
							{
								NovelCTime = CTime[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][e] * pJobSplit[Seq[j]][e];
								pSpeedSelection[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][e] = v;
								break;
							}
						}
						mIdleTime[pStage - 1][mpt] = CTime[pStage - 1][machineAssignment[pStage - 1][i][Seq[j]]][e] * pJobSplit[Seq[j]][e];
					}
				}
			}
			//for the other stages
			for (k = pStage - 2; k >= 0; k--)
			{
				for (int m = 0; m < pMachines[k]; m++)
				{
					for (int j = machineAssignment[k][m].size() - 1; j >= 0; j--)
					{
						for (int e = 0; e < MaxSublotQuantity; e++)
						{
							double NextJobSTime;
							if (j == machineAssignment[k][m].size() - 1)
							{
								NextJobSTime = INT_MAX;
							}
							else
							{
								NextJobSTime = STime[k][machineAssignment[k][m][Seq[j + 1]]][e];
							}

							double NextStageSTime;
							NextStageSTime = TSTime[k + 1][machineAssignment[k][m][Seq[j]]][e];

							int speedLevel = pSpeedSelection[k][machineAssignment[k][m][Seq[j]]][e];

							for (int i = 0; i < speedLevel; i++)
							{
								double difference;
								difference = pUnitTime[k][machineAssignment[k][m][Seq[j]]][i] * pJobSplit[Seq[j]][e] - pUnitTime[k][machineAssignment[k][m][Seq[j]]][speedLevel] * pJobSplit[Seq[j]][e];
								double NovelCTime = CTime[k][machineAssignment[k][m][Seq[j]]][e] * pJobSplit[Seq[j]][e] + difference;
								if (NovelCTime <= min(NextJobSTime, NextStageSTime))
								{
									NovelCTime = CTime[k][machineAssignment[k][m][Seq[j]]][e] * pJobSplit[Seq[j]][e];
									pSpeedSelection[k][machineAssignment[k][m][Seq[j]]][e] = i;
									break;
								}
							}
							mIdleTime[k][mpt] = CTime[k][machineAssignment[k][m][Seq[j]]][e] * pJobSplit[Seq[j]][e];

						}
					}
				}

			}
			//get makespan
			MS = 0;

			for (int i = 0; i < pJob; i++)
			{
				if (CTime[pStage - 1][i][MaxSublotQuantity - 1] > MS)
				{
					MS = CTime[pStage - 1][i][MaxSublotQuantity - 1];
				}

			}
			//get energyconsumption
			TEC = 0.0;
			double EC1 = 0.0;   //energy cost in the processing mode
			double EC2 = 0.0;   //energy cost in the setup mode
			double EC3 = 0.0;   //energy cost in the standby mode

			for (int k = 0; k < pStage; k++)
			{
				for (int i = 0; i < pMachines[k]; i++)
				{
					double totalPT = 0.0;
					for (int j = 0; j < machineAssignment[k][i].size(); j++)
					{
						for (int e = 0; e < MaxSublotQuantity; e++)
						{
							EC1 += pUnitTime[k][machineAssignment[k][i][j]][pSpeedSelection[k][Seq[j]][e]] * pPower[k][machineAssignment[k][i][j]][pSpeedSelection[k][Seq[j]][e]] * pJobSplit[Seq[j]][e];
							totalPT += pUnitTime[k][machineAssignment[k][i][j]][pSpeedSelection[k][Seq[j]][e]];
							EC2 += pSetupTime[k][machineAssignment[k][i][j]] * sp;
							totalPT += pSetupTime[k][machineAssignment[k][i][j]];
						}
					}
					double IdleTime = mIdleTime[k][i] - totalPT;
					EC3 += IdleTime * ip;
				}
			}
			TEC = EC1 + EC2 + EC3;
		}
	}



}




//Individual& Individual::operator=(Individual& obj)
//{
//	this->pJobSeq = obj.pJobSeq;
//	this->pSpeedSelection = obj.pSpeedSelection;
//	this->pAge = obj.pAge;
//	this->MS = obj.MS;
//	this->TEC = obj.TEC;
//
//	this->Fitness = obj.Fitness;
//	this->Fitness = obj.Fitness;
//	this->IndexNeighbor = obj.IndexNeighbor;
//	//this->WeightMS = obj.WeightMS;
//	//this->WeightTEC = obj.WeightTEC;
//	this->CurrentNeighbrohoodIndex = obj.CurrentNeighbrohoodIndex;
//	this->CurrentWeightIndex = obj.CurrentWeightIndex;
//	this->NeighborWeightMS = obj.NeighborWeightMS;
//	this->NeighborWeightTEC = obj.NeighborWeightTEC;
//	return *this;
//}
bool Individual::operator==(Individual& obj)
{
	if (this->MS == obj.MS && this->TEC == obj.TEC)
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool Individual::operator>(Individual& obj)
{
	if (*this == obj)
	{
		return false;
	}
	else
	{
		if ((this->MS <= obj.MS) && (this->TEC <= obj.TEC))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}
bool Individual::operator<(Individual& obj)
{
	if (*this == obj)
	{
		return false;
	}
	else
	{
		if ((this->MS >= obj.MS) && (this->TEC >= obj.TEC))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}
void Individual::getFitnessValue()
{
	//no objective normalization
	//Fitness = 0.0;
	//Fitness += WeightMS * MS;//权重和的方法
	//Fitness += WeightTEC * TEC;
	//objective normalization
	/*Fitness = 0.0;
	Fitness += WeightMS * (MS - MinMS) / (MaxMS - MinMS);
	Fitness += WeightTEC * (TEC - MinTEC) / (MaxTEC - MinTEC);*/
	//Tchebycheff approach objective 切比雪夫方法
	//Tchebycheff approach objective 切比雪夫方法
	double f1;
	f1 = NeighborWeightMS[CurrentWeightIndex] * (MS - MinMS) / (MaxMS - MinMS);
	double f2;
	f2 = NeighborWeightTEC[CurrentWeightIndex] * (TEC - MinTEC) / (MaxTEC - MinTEC);

	if (f1 >= f2)
	{
		Fitness = f1;
	}
	else
	{
		Fitness = f2;
	}
}


Individual::~Individual()
{

}
