#pragma once

#include "Individual.h"
#include "ParetoPoint.h"
#include "CompareIndividual.h"
#include <iomanip>

class NSGAII
{
public:
	NSGAII();
	NSGAII(int _PopulationSize,int _NeighborhoodSize, int _MaximumFuncEvals);
	~NSGAII();
public:
	vector<Individual> ParentPopulation;
	vector<Individual> OffSpringPopulation;
	vector<Individual> ParetoPopulation;
	int PopulationSize;
	int MaximumFuncEvals;
	int MaxinumCPUTime;
	int NeighborhoodSize;
public:
	void InitializePopulation();
	void EvolvePopulation( long time);
	void UpdateParetoPopulation();
	//void IndividualLocalSearch(int k);
	void OutputParetoPopulation(int ins, double costTime);
	//void ConductImprovementStrategy();
	void run(int ins);
	void  TwoPointsXoverAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd);
};

void NSGAII::TwoPointsXoverAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd)
{
	//for the job sequence

	int a, b, m;
	for (int i = 0; i < pJob; i++)
	{
		ResultInd.pJobSeq[i] = pJob;
	}
	//Initialize the pJobSplit
	ResultInd.pJobSplit.resize(pJob);
	for (int i = 0; i < pJob; i++)
	{
		ResultInd.pJobSplit[i].resize(MaxSublotQuantity);
	}

	//选择位置
	vector<int> Position;
	int number = rand() % pJob;
	for (int i = 0; i < number; i++)
	{
		vector<int>::iterator found;
		int position;
		do
		{
			position = rand() % pJob;
			found = find(Position.begin(), Position.end(), position);
		} while (found != Position.end());
		Position.push_back(position);

	}
	vector<int> tempJobs;
	for (int i = 0; i < Position.size(); i++)
	{
		tempJobs.push_back(Ind1.pJobSeq[Position[i]]);
	}

	vector<int> ScheduledPosition;
	for (int i = 0; i < pJob; i++)
	{
		vector<int>::iterator found = find(tempJobs.begin(), tempJobs.end(), Ind2.pJobSeq[i]);
		if (found == tempJobs.end())
		{
			ResultInd.pJobSeq[i] = Ind2.pJobSeq[i];
			ScheduledPosition.push_back(i);
		}
	}
	vector<int> UnscheduledPosition;
	for (int i = 0; i < pJob; i++)
	{
		vector<int>::iterator found = find(ScheduledPosition.begin(), ScheduledPosition.end(), i);
		if (found == ScheduledPosition.end())
		{
			UnscheduledPosition.push_back(i);
		}
	}

	for (int i = 0; i < UnscheduledPosition.size(); i++)
	{
		for (int j = 0; j < pJob; j++)
		{
			vector<int>::iterator found = find(ResultInd.pJobSeq.begin(), ResultInd.pJobSeq.end(), Ind1.pJobSeq[j]);
			if (found == ResultInd.pJobSeq.end())
			{
				ResultInd.pJobSeq[UnscheduledPosition[i]] = Ind1.pJobSeq[j];
				break;
			}
		}
	}
	//for the job split
	for (int i = 0; i < pJob; i++)
	{

		if (rand() * 1.0 / RAND_MAX < 0.416)
		{
			ResultInd.pJobSplit[i] = Ind1.pJobSplit[i];
		}
		else
		{
			ResultInd.pJobSplit[i] = Ind2.pJobSplit[i];
		}
	}
	//conduct the permutation

	if (rand() * 1.0 / RAND_MAX < 0.2)
	{
		int a = rand() % pJob;
		int b = rand() % pJob;
		while (a == b)
		{
			a = rand() % pJob;
			b = rand() % pJob;
		}

		swap(ResultInd.pJobSeq[a], ResultInd.pJobSeq[b]);
	}
	if (rand() * 1.0 / RAND_MAX < 0.2)
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


		int temp1 = ResultInd.pJobSplit[r][a];
		int temp2 = ResultInd.pJobSplit[r][b];

		if (ResultInd.pJobSplit[r][a] > 0 || ResultInd.pJobSplit[r][b] > 0)
		{
			//定义一个步长
			int step = rand() % 3 + 1;
			double r2 = rand() * 1.0 / RAND_MAX;
			if (r2 < 0.5)
			{
				ResultInd.pJobSplit[r][a] += step;
				ResultInd.pJobSplit[r][b] -= step;
			}
			else
			{
				ResultInd.pJobSplit[r][a] -= step;
				ResultInd.pJobSplit[r][b] += step;
			}

			if (ResultInd.pJobSplit[r][0] == 0 && ResultInd.pJobSplit[r][a] < 0 && ResultInd.pJobSplit[r][b] < 0)
			{
				ResultInd.pJobSplit[r][a] = temp1;
				ResultInd.pJobSplit[r][b] = temp2;
			}

		}
	}
	//for the speed selection
	for (int j = 0; j < pStage; j++)
	{

		for (int i = 0; i < pJob; i++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[j][i][e] =1+ rand() % 5;
			}
		}

		//双点交叉操作  
		do
		{
			a = rand() % pJob;
			b = rand() % pJob;
		} while (a == b);

		if (a > b)   //assure that a<b
		{
			m = a;
			a = b;
			b = m;
		}
		for (int i = a; i <= b; i++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[j][i][e] = Ind1.pSpeedSelection[j][i][e];
			}
		}

		for (int i = 0; i < a; i++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[j][i][e] = Ind2.pSpeedSelection[j][i][e];
			}
		}

		for (int i = b + 1; i < pJob; i++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[j][i][e] = Ind2.pSpeedSelection[j][i][e];
			}
		}
	}
	//mutation
	if (rand() * 1.0 / RAND_MAX < 0.1)
	{
		//insertion-based operator for job sequence
		int a, b, e;
		a = rand() % pJob;
		b = rand() % pJob;
		if (b > a)
		{
			int tempValue = ResultInd.pJobSeq[b];
			for (int i = b; i > a; i--)
			{
				ResultInd.pJobSeq[i] = ResultInd.pJobSeq[i - 1];
			}
			ResultInd.pJobSeq[a] = tempValue;
		}
		else
		{
			int tempValue = ResultInd.pJobSeq[b];
			for (int i = b; i < a; i++)
			{
				ResultInd.pJobSeq[i] = ResultInd.pJobSeq[i + 1];
			}
			ResultInd.pJobSeq[a] = tempValue;
		}

		//randomly select a job and swith its processing speed at a randomly picked stage
		a = rand() % pJob;
		b = rand() % pStage;
		e = rand() % MaxSublotQuantity;
		int tempValue;
		do
		{
			tempValue = rand() % 5;
		} while (tempValue == ResultInd.pSpeedSelection[b][a][e]);

		ResultInd.pSpeedSelection[b][a][e] = tempValue;
	}

	ResultInd.getObjectives2(1);
}

NSGAII::NSGAII()
{
}

NSGAII::NSGAII(int _PopulationSize, int _NeighborhoodSize, int _CPU)
{
	this->PopulationSize = _PopulationSize+1;
	this->NeighborhoodSize = _NeighborhoodSize;
	this->MaximumFuncEvals = _CPU*pJob*pStage;
}

void NSGAII::InitializePopulation()
{
	for (int i = 0; i < PopulationSize; i++)
	{
		Individual p(1);
		ParentPopulation.push_back(p);
	}

	OffSpringPopulation.resize(0);//子代解为0
	ParetoPopulation.resize(0);//
}
void NSGAII::EvolvePopulation( long time )
{
	//while (NumberofEvluations < MaximumFuncEvals)
	while((GetTickCount()-time)<MaxinumCPUTime)
	{
		//generate offspring population based on the parentpopulation
		for (int i = 0; i < PopulationSize; i++)
		{
			int a, b;
			do
			{
				a = rand() % PopulationSize;
				b = rand() % PopulationSize;
			} while (a == b);

			Individual tempIndividual(1);
			//NumberofEvluations = NumberofEvluations - 1;
			TwoPointsXoverAndMutation(ParentPopulation[a], ParentPopulation[b], tempIndividual);
			OffSpringPopulation.push_back(tempIndividual);
		}

		for (int i = 0; i < PopulationSize; i++)
		{
			ParentPopulation[i].n = 0;
			ParentPopulation[i].distance = 0;
		}
		//   快速非支配解排序？
		//generate all the nondominated fronts 
		//for parent population
		int maxN = 0;
		for (int i = 0; i < PopulationSize; i++)
		{
			for (int j = 0; j < PopulationSize; j++)
			{
				if (ParentPopulation[i] < ParentPopulation[j])
				{
					ParentPopulation[i].n++;
				}
			}
			for (int j = 0; j < PopulationSize; j++)
			{
				if (ParentPopulation[i] < OffSpringPopulation[j])
				{
					ParentPopulation[i].n++;
				}
			}
			if (ParentPopulation[i].n > maxN)
			{
				maxN = ParentPopulation[i].n;
			}
		}
		//for offspring population
		for (int i = 0; i < PopulationSize; i++)
		{
			for (int j = 0; j < PopulationSize; j++)
			{
				if (OffSpringPopulation[i] < ParentPopulation[j])
				{
					OffSpringPopulation[i].n++;
				}
			}
			for (int j = 0; j < PopulationSize; j++)
			{
				if (OffSpringPopulation[i] < OffSpringPopulation[j])
				{
					OffSpringPopulation[i].n++;
				}
			}
			if (OffSpringPopulation[i].n > maxN)
			{
				maxN = OffSpringPopulation[i].n;
			}
		}
		//generate nondominated fronts method(1)
		vector<vector<Individual>> p;
		p.resize(maxN + 1);
		for (int i = 0; i < PopulationSize; i++)
		{
			p[ParentPopulation[i].n].push_back(ParentPopulation[i]);
		}
		for (int i = 0; i < PopulationSize; i++)
		{
			p[OffSpringPopulation[i].n].push_back(OffSpringPopulation[i]);
		}
		//calculate the crowding distance in P[i]计算拥挤距离
		for (int i = 0; i < p.size(); i++)
		{
			if (p[i].size() > 0)
			{
				//compare MS
				sort(p[i].begin(), p[i].end(), CompareMS());//目标函数进行排序
				//给边界点之外的点赋值
				for (int j = 1; j < p[i].size() - 1; j++)
				{
					p[i][j].distance += abs(p[i][j + 1].MS - p[i][j - 1].MS);//相邻解绝对值差
				}
				//给边界点赋值 赋值大值保证能被选中
				p[i][0].distance += 10000;
				p[i][p[i].size() - 1].distance += 10000;
				//compare TEC
				sort(p[i].begin(), p[i].end(), CompareTEC());
				//给边界点之外的点赋值
				for (int j = 1; j < p[i].size() - 1; j++)
				{
					p[i][j].distance += abs(p[i][j + 1].TEC - p[i][j - 1].TEC);
				}
				//给边界点赋值 赋值大值保证能被选中
				p[i][0].distance += 10000;
				p[i][p[i].size() - 1].distance += 10000;
				sort(p[i].begin(), p[i].end(), CompareDistance());
			}
		}
		//update the parent population
		ParentPopulation.clear();//父代解
		OffSpringPopulation.clear();//子代解
		bool isFull = false;
		for (int i = 0; i < p.size(); i++)
		{
			if (!isFull)
			{
				for (int j = 0; j < p[i].size(); j++)
				{
					ParentPopulation.push_back(p[i][j]);
					if (ParentPopulation.size() >= PopulationSize)
					{
						isFull = true;
						break;
					}
				}
			}
		}
		//update ParetoPopulation
		UpdateParetoPopulation();
	}
}

void NSGAII::UpdateParetoPopulation()
{
	//update the pareto population
	for (int i = 0; i < PopulationSize; i++)
	{
		//remove from ParetoPopulation dominated by this
		for (int n = 0; n < ParetoPopulation.size(); n++)
		{
			if (ParentPopulation[i] > ParetoPopulation[n])
			{
				ParetoPopulation.erase(ParetoPopulation.begin() + n);
				n--;
			}
		}
		//Add this to ParetoPopulation if no points dominate this
		bool  isDominated = false;
		for (vector<Individual>::iterator it = ParetoPopulation.begin(); it != ParetoPopulation.end(); it++)
		{
			if (*it > ParentPopulation[i] || *it == ParentPopulation[i])
			{
				isDominated = true;
				break;
			}
		}

		if (!isDominated)
		{
			ParetoPopulation.push_back(ParentPopulation[i]);
		}
	}
}
void NSGAII::OutputParetoPopulation(int ins, double costTime)
{
	/*for (int i = 0; i < ParetoPopulation.size(); i++)
	{
	cout << ParetoPopulation[i].MS << " " << ParetoPopulation[i].MF << " " << ParetoPopulation[i].TF << endl;
	}*/

	char outfile[40];
	sprintf_s(outfile, "results.nsga2\\%d_%d_%d_%d.txt", pJob, pStage, pType, ins);
	fstream fout(outfile, fstream::out | fstream::app);
	fout << ParetoPopulation.size() << " " << costTime << endl;
	/*fout << "MS=[";*/
	for (int i = 0; i < ParetoPopulation.size(); i++)
	{
		fout << ParetoPopulation[i].MS << " ";
	}
	fout << endl;
	/*fout << "];" << endl;*/


	/*fout << "TF=[";*/
	for (int i = 0; i < ParetoPopulation.size(); i++)
	{
		fout << ParetoPopulation[i].TEC << " ";
	}
	fout << endl;
	/*fout << "];" << endl;*/

	fout.close();
}

void NSGAII::run(int ins)
{
	long initTime, finalTime;
	double costTime;
	initTime = GetTickCount();
	InitializePopulation();
	UpdateParetoPopulation();
	EvolvePopulation(initTime);
	//ConductImprovementStrategy();
	finalTime = GetTickCount();
	costTime = (finalTime - initTime) / 1000.0;
	OutputParetoPopulation(ins, costTime);
}


NSGAII::~NSGAII()
{
}



