#pragma once
#include "Individual.h"
#include "ParetoPoint.h"
#include "CompareIndividual.h"
#include <iomanip>
class MOCGWO
{
public:
	MOCGWO();
	MOCGWO(int _PopulationSize, int _CPU);
	~MOCGWO();
public:
	void InitializePopulation();
	void InitializeGridsAndFitness();
	void EvolvePopulation(long time);
	void UpdateParetoPopulation();
	void UpdateParetoPopulation(Individual& Ind);
	void IndividualLocalSearch(int k);
	//void OutputParetoPopulation(int ins, double costTime);
	//void ConductImprovementStrategy();
	void run(int ins);
	//void getFitnessForIndividual(Individual& Ind);

	//void static TwoPointsXoverAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd);


public:
	vector<Individual>  Population;   //current population
	vector<Individual>  ParetoPopulation;   //pareto points
public:
	int PopulationSize;    //the current population size
	int NeighborhoodSize;  //the neighborhood size
	int MaximumCPUTime;  //the maximum function evaluations
	int AgeTrail;          //age trial
	int VNSTrial;         //VNS trial
	bool IsVD;            //是否采用VD策略
	int RestartStrategy;   //重启策略
public:
	//for the weight vectors of two objectives
	int H;
	vector<double> W1;
	vector<double> W2;

public:
	void TwoPointsXoverAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd);
	//void UpdateReferencePoint(Individual &Ind);
	void OutputParetoPopulation(int ins, double costTime);

public:
	ParetoPoint RP;  //reference point 
};

MOCGWO::~MOCGWO()
{

}
MOCGWO::MOCGWO(int _PopulationSize, int _CPU)
{
	this->PopulationSize = _PopulationSize;
	this->MaximumCPUTime = _CPU * pJob * pStage;
}
//void MOCGWO::InitializeGridsAndFitness()
//{
//	for (int i = 0; i < PopulationSize; i++)
//	{
//		vector<int>::iterator it;
//		for (int j = 0; j < 10; j++)
//		{
//			int tempIndex;
//			do
//			{
//				tempIndex = rand() % PopulationSize;
//				it = find(Population[i].IndexNeighbor.begin(), Population[i].IndexNeighbor.end(), tempIndex);
//			} while (it != Population[i].IndexNeighbor.end());
//
//			Population[i].IndexNeighbor.push_back(tempIndex);
//		}
//	}
//}


void MOCGWO::InitializePopulation()
{
	//basic information
	for (int i = 0; i < PopulationSize; i++)
	{
		Individual p(1);
		Population.push_back(p);
	}
	////add weight information
	//for (int i = 0; i < PopulationSize; i++)
	//{
	//	Population[i].WeightMS = W1[i];
	//	Population[i].WeightTEC = W2[i];
	//}
	//get the fitness informaiton
	for (int i = 0; i < PopulationSize; i++)
	{
		Population[i].getFitnessValue();
	}
	////NondominatedSorting();
	////NondominatedSorting();

	//update the ParetoPopulation
	for (int i = 0; i < PopulationSize; i++)
	{
		UpdateParetoPopulation(Population[i]);
	}
	//initialize the indexNeighbor and neighborWeight
	for (int i = 0; i < PopulationSize; i++)
	{
		//Pair<double>* ch = new Pair<double>[PopulationSize];
		Pair1<double>* ch = new Pair1<double>[PopulationSize];
		for (int j = 0; j < PopulationSize; j++)
		{
			ch[j].dim = j;
			double tempValue = 0;
			tempValue += (Population[i].WeightMS - Population[j].WeightMS) * (Population[i].WeightMS - Population[j].WeightMS);
			tempValue += (Population[i].WeightTEC - Population[j].WeightTEC) * (Population[i].WeightTEC - Population[j].WeightTEC);

			ch[j].weightMS = Population[j].WeightMS;
			ch[j].weightTEC = Population[j].WeightTEC;
			ch[j].value = sqrt(tempValue);
		}
		sort(ch, ch + PopulationSize, PairLess1<double>());  //从小到大打排序a
		Population[i].NeighborWeightMS.clear();
		Population[i].NeighborWeightTEC.clear();
		for (int k = 0; k < NeighborhoodSize; k++)
		{
			//Population[i].IndexNeighbor.push_back(ch[k].dim);
			Population[i].NeighborWeightMS.push_back(ch[k].weightMS);
			Population[i].NeighborWeightTEC.push_back(ch[k].weightTEC);
		}
		Population[i].IndexNeighbor.clear();

		for (int k = 0; k < NeighborhoodSize; k++)
		{
			Population[i].IndexNeighbor.push_back(ch[k].dim);
		}

		delete[] ch;
	}
}
void MOCGWO::EvolvePopulation(long time)
{
	//int generations = 0;
	long InitTime = GetTickCount();
	//while ((GetTickCount() - time) < MaximumCPUTime)
	while ((GetTickCount() - time) < MaximumCPUTime)
	{
		//generations++;
		for (int i = 0; i < PopulationSize; i++)
		{
			//select three best solutions from the neighbors according to the the fitness 
			int FirstIndex;
			int SecondIndex;
			int ThirdIndex;
			//select the first Index
			FirstIndex = 0;
			for (int j = 1; j < 10; j++)
			{
				if (Population[j].Fitness < Population[FirstIndex].Fitness)

					FirstIndex = j;
			}

			//select the second Index
			do
			{
				SecondIndex = rand() % 10;
			} while (SecondIndex == FirstIndex);

			for (int j = 0; j < 10; j++)
			{
				if (j != FirstIndex)
				{
					if (Population[j].Fitness < Population[SecondIndex].Fitness)
					{
						SecondIndex = j;
					}
				}
			}
			//select the third Index
			do
			{
				ThirdIndex = rand() % 10;
			} while (ThirdIndex == FirstIndex || ThirdIndex == SecondIndex);

			for (int j = 0; j < 10; j++)
			{
				if (j != FirstIndex && j != SecondIndex)
				{
					if (Population[j].Fitness < Population[ThirdIndex].Fitness)
					{
						ThirdIndex = j;
					}
				}
			}
			//conduct the search operator
			//select the collaborative index
			int mapIndex;
			double r = rand() * 1.0 / RAND_MAX;
			if (r < 1.0 / 3)
			{
				mapIndex = FirstIndex;
			}
			else if (r >= 1.0 / 3 && r < 2.0 / 3)
			{
				mapIndex = SecondIndex;
			}
			else
			{
				mapIndex = ThirdIndex;
			}

			Individual tempIndividual;
			tempIndividual = Population[i];

			TwoPointsXoverAndMutation(Population[i], Population[mapIndex], tempIndividual);

			UpdateParetoPopulation(tempIndividual);

			if (tempIndividual > Population[i])
			{
				Population[i] = tempIndividual;
			}

			//conduct the VNS procedure
			for (int i = 0; i < 4; i++)
			{
				//select a solution from the population
				int tempIndex = rand() % PopulationSize;

				//congduct the varaible neighborhood search
				int k = 1;
				while (k <= 8)
				{
					Individual Neighbor;
					Neighbor = Population[tempIndex];
					Population[tempIndex].getNeigbors(Neighbor, k);

					UpdateParetoPopulation(Neighbor);
					//get the fitness value for the Neighbor
					 //getFitnessForIndividual(Neighbor);

					if (Neighbor > Population[tempIndex])
					{
						k = 1;
						Population[tempIndex] = Neighbor;
					}
					else if (Neighbor < Population[tempIndex])
					{
						k++;
					}

					if (Neighbor.Fitness < Population[tempIndex].Fitness)
					{
						k = 1;
						Population[tempIndex] = Neighbor;
					}
					else
					{
						k++;
					}
				}
			}
			//Update the fitness value for the popualtion
				//Calcuate the density value
			for (int i = 0; i < PopulationSize; i++)
			{
				//get the nearest neighbor with the Euclidean value
				double minEuclidean = INT_MAX;
				for (int j = 0; j < 10; j++)
				{
					double tempValue = sqrt((Population[i].MS - Population[Population[i].IndexNeighbor[j]].MS) * (Population[i].MS - Population[Population[i].IndexNeighbor[j]].MS)
						+ (Population[i].TEC - Population[Population[i].IndexNeighbor[j]].TEC) * (Population[i].TEC - Population[Population[i].IndexNeighbor[j]].TEC));
					if (tempValue < minEuclidean)
					{
						minEuclidean = tempValue;
					}
				}
				//Population[i].DensityValue = 1.0 / (minEuclidean + 2);
			}

			////Initialize the raw fitness value
			//for (int i = 0; i < PopulationSize; i++)
			//{
			//	for (int j = 0; j < PopulationSize; j++)
			//	{
			//		if (Population[i] > Population[j])
			//		{
			//			Population[i].StrengthValue++;
			//		}
			//	}
			//}

			//for (int i = 0; i < PopulationSize; i++)
			//{
			//	for (int j = 0; j < PopulationSize; j++)
			//	{
			//		if (Population[j] > Population[i])
			//		{
			//			Population[i].RawValue += Population[j].StrengthValue;
			//		}
			//	}
			//}

			//for (int i = 0; i < PopulationSize; i++)
			//{
			//	Population[i].Fitness = Population[i].DensityValue + Population[i].RawValue;
			//}


		}
	}
}



void MOCGWO::UpdateParetoPopulation()
{
	//update the pareto population
	for (int i = 0; i < PopulationSize; i++)
	{
		//remove from ParetoPopulation dominated by this
		for (int n = 0; n < ParetoPopulation.size(); n++)
		{
			if (Population[i] > ParetoPopulation[n])
			{
				ParetoPopulation.erase(ParetoPopulation.begin() + n);
				n--;
			}
		}
		//Add this to ParetoPopulation if no points dominate this
		bool  isDominated = false;
		for (vector<Individual>::iterator it = ParetoPopulation.begin(); it != ParetoPopulation.end(); it++)
		{
			if (*it > Population[i] || *it == Population[i])
			{
				isDominated = true;
				break;
			}
		}

		if (!isDominated)
		{
			ParetoPopulation.push_back(Population[i]);
		}
	}
}
void MOCGWO::UpdateParetoPopulation(Individual& Ind)
{
	//remove from ParetoPopulation dominated by this
	for (int n = 0; n < ParetoPopulation.size(); n++)
	{
		if (Ind > ParetoPopulation[n])
		{
			ParetoPopulation.erase(ParetoPopulation.begin() + n);
			n--;
		}
	}
	//Add this to ParetoPopulation if no points dominate this
	bool  isDominated = false;
	for (vector<Individual>::iterator it = ParetoPopulation.begin(); it != ParetoPopulation.end(); it++)
	{
		if (*it > Ind || *it == Ind)
		{
			isDominated = true;
			break;
		}
	}

	if (!isDominated)
	{
		ParetoPopulation.push_back(Ind);
	}
}
void  MOCGWO::TwoPointsXoverAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd)
{
	//for the job sequence
		//双点交叉操作
	ResultInd.pJobSeq.resize(pJob);
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
	ResultInd.pSpeedSelection.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		ResultInd.pSpeedSelection[k].resize(pJob);
	}
	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[k][j].resize(MaxSublotQuantity);
			}
		}
	}
	//for the speed selection
	for (int k = 0; k < pStage; k++)
	{
		int a, b, m;
		int temp;

		for (int j = 0; j < pJob; j++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[k][j][e] = 5;
			}
		}
		//做双点交叉操作
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
				ResultInd.pSpeedSelection[k][i][e] = Ind1.pSpeedSelection[k][i][e];
			}
		}

		for (int i = 0; i < a; i++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[k][i][e] = Ind2.pSpeedSelection[k][i][e];
			}
		}

		for (int i = b + 1; i < pJob; i++)
		{
			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				ResultInd.pSpeedSelection[k][i][e] = Ind2.pSpeedSelection[k][i][e];
			}
		}

	}
	//conduct mutation
	if (rand() * 1.0 / RAND_MAX < 0.1)
	{
		int a, b;
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
		int e = rand() % MaxSublotQuantity;
		int tempValue;
		do
		{
			tempValue = rand() % 5;
		} while (tempValue == ResultInd.pSpeedSelection[b][a][e]);

		ResultInd.pSpeedSelection[b][a][e] = tempValue;

	}
	ResultInd.getObjectives2(1);
}
void MOCGWO::run(int ins)
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
void MOCGWO::OutputParetoPopulation(int ins, double costTime)
{
	sort(ParetoPopulation.begin(), ParetoPopulation.end(), CompareMS());
	/*for (int i = 0; i < ParetoPopulation.size(); i++)
	{
	cout << ParetoPopulation[i].MS << " " << ParetoPopulation[i].MF << " " << ParetoPopulation[i].TF << endl;
	}*/

	char outfile[40];
	sprintf_s(outfile, "results.MOCGWO\\%d_%d_%d_%d.txt", pJob, pStage, pType, ins);
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
