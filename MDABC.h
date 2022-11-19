#pragma once
#include "Individual.h"
#include "ParetoPoint.h"
#include "CompareIndividual.h"
#include <iomanip>


class MDABC
{
public:
	MDABC();
	MDABC(int _H, int _NeighborhoodSize,  int _NeighborhoodTrial, int _AgeTrial, bool _IsVD,  int _MaximumCPU);
	~MDABC();
public:
	void Gen2WV_NBI(); //generate weight vectors for two objectives
	void UpdateParetoPouplation(Individual& Ind);
	void InitializePopulation();  //initialize the population
	void EvolvePopulation(long time);      //Evolve the population
	void NondominatedSorting();     //non-dominated sorting
	void ConductImprovementStrategy();

	void run(int ins);


public:
	vector<Individual>  Population;   //current population
	vector<Individual>  ParetoPopulation;   //pareto points
	vector<unsigned short>bestjodseq;
public:
	int PopulationSize;    //the current population size
	int NeighborhoodSize;  //the neighborhood size
	int NeighborhoodTrial;   //neighborhood trial
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
	void PBXAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd);
	//void UpdateReferencePoint(Individual &Ind);
	void OutputParetoPopulation(int ins, double costTime);

public:
	ParetoPoint RP;  //reference point 
};


MDABC::MDABC()
{

}

MDABC::~MDABC()
{

}

MDABC::MDABC(int _H, int _NeighborhoodSize,  int _NeighborhoodTrial, int _AgeTrial, bool _IsVD,  int _MaximumCPU)
{
	this->H = _H;
	this->NeighborhoodSize = _NeighborhoodSize;
	this->NeighborhoodTrial = _NeighborhoodTrial;
	this->AgeTrail = _AgeTrial;
	//this->VNSTrial = _VNSTrial;
	this->IsVD = _IsVD;
	
	this->MaximumCPUTime = _MaximumCPU * pJob * pStage;

	this->PopulationSize = H + 1; //149->150  199->200  249->250
}
void MDABC::Gen2WV_NBI()
{
	W1.resize(PopulationSize);
	W2.resize(PopulationSize);
	for (int i = 0, j = H; i <= H, j >= 0; i++, j--)
	{
		W1[i] = i * 1.0 / H;
		W2[i] = j * 1.0 / H;
	}
	//int a;
	//调试用
	//for (int i = 0; i < PopulationSize; i++)
	//{
	//	cout << W1[i] + W2[i] << endl;
	//}

}

void MDABC::run(int ins)
{
	long initTime, finalTime;
	double costTime;
	initTime = GetTickCount();
	Gen2WV_NBI();
	InitializePopulation();
	EvolvePopulation(initTime);
	//ConductImprovementStrategy();
	finalTime = GetTickCount();
	costTime = (finalTime - initTime) / 1000.0;
	OutputParetoPopulation(ins, costTime);
}

void MDABC::InitializePopulation()
{
	//basic information
	for (int i = 0; i < PopulationSize; i++)
	{
		Individual p(1);
		Population.push_back(p);
	}

	//add weight information
	for (int i = 0; i < PopulationSize; i++)
	{
		Population[i].WeightMS = W1[i];
		Population[i].WeightTEC = W2[i];
	}
	//get the fitness informaiton
	for (int i = 0; i < PopulationSize; i++)
	{
		Population[i].getFitnessValue();
	}
	//NondominatedSorting();

	////update the reference point
	//RP.MS = INT_MAX * 1.0;
	//RP.TEC = INT_MAX * 1.0;
	//for (int i = 0; i < PopulationSize; i++)
	//{
	//	UpdateReferencePoint(Population[i]);
	//}
	//update the ParetoPopulation更新pa雷托
	for (int i = 0; i < PopulationSize; i++)
	{
		UpdateParetoPouplation(Population[i]);
	}
	//initialize the indexNeighbor
	for (int i = 0; i < PopulationSize; i++)
	{
		Pair1<double>* ch = new Pair1<double>[PopulationSize];
		for (int j = 0; j < PopulationSize; j++)
		{
			ch[j].dim = j;
			double tempValue = 0;
			tempValue += (Population[i].WeightMS - Population[j].WeightMS) * (Population[i].WeightMS - Population[j].WeightMS);
			tempValue += (Population[i].WeightTEC - Population[j].WeightTEC) * (Population[i].WeightTEC - Population[j].WeightTEC);
			ch[j].weightMS = Population[j].WeightMS;
			ch[j].weightTEC = Population[j].WeightTEC;
			

			if (IsVD)
			{
				ch[j].value = sqrt(tempValue) * (Population[j].n + 1);
			}
			else
			{
				ch[j].value = sqrt(tempValue);
			}
		}
		sort(ch, ch + PopulationSize, PairLess1<double>());  //从小到大打排xu
		//Population[i].IndexNeighbor.clear();
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
//void MDABC::UpdateReferencePoint(Individual &Ind)
//{
//	if (Ind.MS < RP.MS)
//	{
//		RP.MS = Ind.MS;
//	}
//	if (Ind.TEC < RP.TEC)
//	{
//		RP.TEC = Ind.TEC;
//	}
//}

void MDABC::ConductImprovementStrategy()
{
	/*for (int i = 0; i < ParetoPopulation.size(); i++)
	{
		ParetoPopulation[i].getObjectives4(1);
	}*/

	vector<Individual> tempPopulation = ParetoPopulation;

	ParetoPopulation.clear();

	//update the pareto population
	for (int i = 0; i < tempPopulation.size(); i++)
	{
		//remove from ParetoPopulation dominated by this
		for (int n = 0; n < ParetoPopulation.size(); n++)
		{
			if (tempPopulation[i] > ParetoPopulation[n])
			{
				ParetoPopulation.erase(ParetoPopulation.begin() + n);
				n--;
			}
		}
		//Add this to ParetoPopulation if no points dominate this
		bool  isDominated = false;
		for (vector<Individual>::iterator it = ParetoPopulation.begin(); it != ParetoPopulation.end(); it++)
		{
			if (*it > tempPopulation[i] || *it == tempPopulation[i])
			{
				isDominated = true;
				break;
			}
		}

		if (!isDominated)
		{
			ParetoPopulation.push_back(tempPopulation[i]);
		}
	}
}

void MDABC::PBXAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd)
{
	ResultInd.pJobSeq.resize(pJob);
	//for the job sequence

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
	//initialize the pspeedselection
	ResultInd.pSpeedSelection.resize(pStage);
	for (int k = 0; k < pStage; k++)
	{
		ResultInd.pSpeedSelection[k].resize(pJob);
	}
	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{
			ResultInd.pSpeedSelection[k][j].resize(MaxSublotQuantity);
		}
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


	for (int i = 0; i < Position.size(); i++)
	{
		ResultInd.pJobSeq[Position[i]] = Ind1.pJobSeq[Position[i]];
	}

	vector<int> UnscheduledPosition;
	for (int i = 0; i < pJob; i++)
	{
		vector<int>::iterator found = find(Position.begin(), Position.end(), i);
		if (found == Position.end())
		{
			UnscheduledPosition.push_back(i);
		}
	}

	for (int i = 0; i < UnscheduledPosition.size(); i++)
	{
		for (int j = 0; j < pJob; j++)
		{
			vector<int>::iterator found = find(ResultInd.pJobSeq.begin(), ResultInd.pJobSeq.end(), Ind2.pJobSeq[j]);
			if (found == ResultInd.pJobSeq.end())
			{
				ResultInd.pJobSeq[UnscheduledPosition[i]] = Ind2.pJobSeq[j];
				break;
			}
		}

	}

	//for the job split
	for (int i = 0; i < pJob; i++)
	{
			if (rand() * 1.0 / RAND_MAX < 0.5)
			{
				ResultInd.pJobSplit[i] = Ind1.pJobSplit[i];
			}
			else
			{
				ResultInd.pJobSplit[i] = Ind2.pJobSplit[i];
			}
		
	}
	//for the speed selection
	for (int k = 0; k < pStage; k++)
	{
		for (int j = 0; j < pJob; j++)
		{

			for (int e = 0; e < MaxSublotQuantity; e++)
			{
				if (rand() * 1.0 / RAND_MAX < 0.5)
				{
					ResultInd.pSpeedSelection[k][j][e] = Ind1.pSpeedSelection[k][j][e];
				}
				else
				{
					ResultInd.pSpeedSelection[k][j][e] = Ind2.pSpeedSelection[k][j][e];
				}
			}
		}
	//	//做双点交叉操作
	//	do
	//	{
	//		a = rand() % pJob;
	//		b = rand() % pJob;
	//	} while (a == b);

	//	if (a > b)   //assure that a<b
	//	{
	//		m = a;
	//		a = b;
	//		b = m;
	//	}
	//	for (int i = a; i <= b; i++)
	//	{
	//		for (int e = 0; e < MaxSublotQuantity; e++)
	//		{
	//			ResultInd.pSpeedSelection[k][i][e] = Ind1.pSpeedSelection[k][i][e];
	//		}
	//	}

	//	for (int i = 0; i < a; i++)
	//	{
	//		for (int e = 0; e < MaxSublotQuantity; e++)
	//		{
	//			ResultInd.pSpeedSelection[k][i][e] = Ind2.pSpeedSelection[k][i][e];
	//		}
	//	}

	//	for (int i = b + 1; i < pJob; i++)
	//	{
	//		for (int e = 0; e < MaxSublotQuantity; e++)
	//		{
	//			ResultInd.pSpeedSelection[k][i][e] = Ind2.pSpeedSelection[k][i][e];
	//		}
	//	}

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

		//randomly select a job and swith its processing speed at a randomly picked stage
		/*if (rand() * 1.0 / RAND_MAX < 0.2)
		{
			 int a = rand() % pJob;
			 int b = rand() % pStage;
			 int e = rand() % MaxSublotQuantity;
			int tempValue;
			do
			{
				tempValue = rand() % 5;
			} while (tempValue == ResultInd.pSpeedSelection[b][a][e]);

			ResultInd.pSpeedSelection[b][a][e] = tempValue;
		}*/
		for (int k = 0;k < pStage;k++)
		{
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


				int temp1 = ResultInd.pSpeedSelection[k][r][a];
				int temp2 = ResultInd.pSpeedSelection[k][r][b];

				if (ResultInd.pSpeedSelection[k][r][a] > 0 || ResultInd.pSpeedSelection[k][r][b] > 0)
				{
					//定义一个步长
					int step = rand() % 3 + 1;
					double r2 = rand() * 1.0 / RAND_MAX;
					if (r2 < 0.5)
					{
						ResultInd.pSpeedSelection[k][r][a] += step;
						ResultInd.pSpeedSelection[k][r][b] -= step;
					}
					else
					{
						ResultInd.pSpeedSelection[k][r][a] -= step;
						ResultInd.pSpeedSelection[k][r][b] += step;
					}

					if (ResultInd.pSpeedSelection[k][r][0] == 0 )
					{
						ResultInd.pSpeedSelection[k][r][a] = temp1;
						ResultInd.pSpeedSelection[k][r][b] = temp2;
					}

				}
			}
		}
	}
	   ResultInd.getObjectives2(1);
}

void MDABC::EvolvePopulation(long time)
{
	//int generations = 0;
	long InitTime = GetTickCount();
	while((GetTickCount() - time) < MaximumCPUTime)
	{//while (generations < 1000)
	
		//generations++;
		//employed bee phase雇佣蜂
		for (int i = 0; i < PopulationSize; i++)
		{
			//Individual *neighbor = Population[i].getNeigbors(Population[i].NS);
			Individual Neighbor;
			Neighbor = Population[i];
			Population[i].getNeigbors(Neighbor, Population[i].NS);
			UpdateParetoPouplation(Neighbor);
			//UpdateReferencePoint(*neighbor);
			if (Neighbor.Fitness < Population[i].Fitness)
			{
				Population[i] = Neighbor;
				Population[i].NS = 1;
				Population[i].getFitnessValue();
				Population[i].pAge = 0;

			}
			else
			{
				Population[i].pAge++;
				Population[i].NS++;
			}
			//delete Neighbor;

			if (Population[i].NS > 8)
			{

				Population[i].NS = 1;
			}
		}
		//only re-initialzie the n value
		//NondominatedSorting();
		//re-initialize the indexNeighbor
		for (int i = 0; i < PopulationSize; i++)
		{
			Pair1<double>* ch = new Pair1<double>[PopulationSize];
			for (int j = 0; j < PopulationSize; j++)
			{
				ch[j].dim = j;
				double tempValue = 0;
				tempValue += (Population[i].WeightMS - Population[j].WeightMS) * (Population[i].WeightMS - Population[j].WeightMS);
				tempValue += (Population[i].WeightTEC - Population[j].WeightTEC) * (Population[i].WeightTEC - Population[j].WeightTEC);
				ch[j].weightMS = Population[j].WeightMS;
				ch[j].weightTEC = Population[j].WeightTEC;
				if (IsVD)
				{
					ch[j].value = sqrt(tempValue) * (Population[j].n + 1);
				}
				else
				{
					ch[j].value = sqrt(tempValue);
				}

			}

			sort(ch, ch + PopulationSize, PairLess1<double>());  //从小到大打排序

			Population[i].IndexNeighbor.clear();


			for (int k = 0; k < NeighborhoodSize; k++)
			{
				Population[i].IndexNeighbor.push_back(ch[k].dim);
			}

			delete[] ch;
		}
		//onlooer bee 
		for (int i = 0; i < PopulationSize; i++)
		{
			int r[3];
			do
			{
				r[0] = rand() % PopulationSize;
				r[1] = rand() % PopulationSize;
				r[2] = rand() % PopulationSize;
			} while (r[0] == r[1] || r[0] == r[2] || r[1] == r[2]);

			double topsis_pdistance[3];      //topsis_distance
			double topsis_ndistance[3];      //topsis_ndistance
			double topsis_rdistance[3];       //relative distance

			for (int j = 0; j < 3; j++)
			{
				topsis_ndistance[j] = sqrt((Population[r[j]].MS - MaxMS) * (Population[r[j]].MS - MaxMS) + (Population[r[j]].TEC - MaxTEC) * (Population[r[j]].TEC - MaxTEC));
				topsis_pdistance[j] = sqrt((Population[r[j]].MS - MinMS) * (Population[r[j]].MS - MinMS) + (Population[r[j]].TEC - MinTEC) * (Population[r[j]].TEC - MinTEC));
			}

			for (int j = 0; j < 3; j++)
			{
				topsis_rdistance[j] = topsis_ndistance[j] / (topsis_ndistance[j] + topsis_pdistance[j]);
			}

			int maxIndex = r[0];
			for (int j = 1; j < 3; j++)
			{
				if (topsis_rdistance[j] > topsis_rdistance[maxIndex])
				{
					maxIndex = r[j];
				}
			}
			Individual tempInd = Population[maxIndex];//父代解
			tempInd.pAge = 0;
			int r3 = rand() % NeighborhoodSize;//随机选邻域解
			int a = Population[maxIndex].IndexNeighbor[r3];
			PBXAndMutation(Population[maxIndex], Population[a], tempInd);//
			//update the reference point
			//UpdateReferencePoint(tempInd);
			//update the ParetoPopulation
			UpdateParetoPouplation(tempInd);
			int select[2];
			select[0] = maxIndex;
			select[1] = a;
			//update the selected solution with the global replacement strategy
			for (int k = 0; k < 2; k++)
			{
				tempInd.NeighborWeightMS = Population[select[k]].NeighborWeightMS;
				tempInd.NeighborWeightTEC = Population[select[k]].NeighborWeightTEC;

				if (Population[select[k]].CurrentWeightIndex == 0)
				{
					tempInd.CurrentWeightIndex = 0;
					tempInd.getFitnessValue();

					if (tempInd.Fitness < Population[select[k]].Fitness)
					{
						Population[select[k]].pJobSeq = tempInd.pJobSeq;
						Population[select[k]].pJobSplit = tempInd.pJobSplit;
						Population[select[k]].pSpeedSelection = tempInd.pSpeedSelection;
						Population[select[k]].MS = tempInd.MS;
						Population[select[k]].TEC = tempInd.TEC;
						Population[select[k]].Fitness = tempInd.Fitness;
						Population[select[k]].NS = 1;
						Population[select[k]].NSCF = 0;
						Population[select[k]].pAge = 0;
						Population[select[k]].neighborhoodAge = 0;
					}
					else
					{
						Population[select[k]].pAge++;
					}
				}
				else
				{
					int tempWeight = Population[select[k]].CurrentWeightIndex;
					double tempFitness = Population[select[k]].Fitness;
					bool IsImprovement = false;

					for (int g = 0; g < NeighborhoodSize; g++)
					{
						Population[select[k]].CurrentWeightIndex = g;
						Population[select[k]].getFitnessValue();
						tempInd.CurrentWeightIndex = g;
						tempInd.getFitnessValue();
						if (tempInd.Fitness < Population[select[k]].Fitness)
						{
							Population[select[k]].pJobSeq = tempInd.pJobSeq;
							Population[select[k]].pJobSplit= tempInd.pJobSplit;
							Population[select[k]].pSpeedSelection = tempInd.pSpeedSelection;
							Population[select[k]].MS = tempInd.MS;
							Population[select[k]].TEC = tempInd.TEC;
							Population[select[k]].CurrentWeightIndex = 0;
							Population[select[k]].getFitnessValue();
							Population[select[k]].NS = 1;
							Population[select[k]].NSCF = 0;

							IsImprovement = true;
							break;
						}
					}

					if (!IsImprovement)
					{
						Population[select[k]].CurrentWeightIndex = tempWeight;
						Population[select[k]].Fitness = tempFitness;
						Population[select[k]].pAge++;
					}
					else
					{
						Population[select[k]].pAge = 0;
					}
				}

			}

		}






		//	//update the selected solutions
		//	tempInd.getFitnessValue();
		//	if (tempInd.Fitness < Population[maxIndex].Fitness)
		//	{
		//		Population[maxIndex] = tempInd;
		//		Population[maxIndex].NS = 1;
		//		Population[maxIndex].NSCF = 0;
		//	}
		//	//update the neighboring solution
		//	tempInd.WeightMS = Population[a].WeightMS;
		//	tempInd.WeightTEC = Population[a].WeightTEC;
		//	tempInd.getFitnessValue();



		//	if (tempInd.Fitness < Population[a].Fitness)
		//	{
		//		Population[a] = tempInd;
		//		Population[a].NS = 1;
		//		Population[a].NSCF = 0;
		//	}
		//}

		//Scout bees
		/*for (int i = 0; i < PopulationSize; i++)
		{
			if (Population[i].pAge > AgeTrail)
			{
				Population[i].CurrentWeightIndex++;
				if (Population[i].CurrentWeightIndex < NeighborhoodSize)
				{
					Population[i].getFitnessValue();
					Population[i].getObjectives2(1);
					Population[i].pAge = 0;
				}
			}
			else
			{
				Population[i].CurrentWeightIndex = 0;
				Population[i].pAge = 0;
				Population[i].conductPermution();
				Population[i].getFitnessValue();
				Population[i].getObjectives2(1);
			}



		}
		 */
		 
		for (int i = 0; i < PopulationSize; i++)
		{
			//if (RestartStrategy == 1)
			//{
			//	if (Population[i].pAge > AgeTrail)
			//	{
			//		Population[i].CurrentWeightIndex++;
			//		bool IsChanged = false;
			//		for (int j = 0; j < NeighborhoodSize; j++)
			//		{
			//			if (Population[i].CurrentWeightIndex < NeighborhoodSize)
			//			{
			//				IsChanged = true;

			//				//交换他们的解
			//				vector<int> tempJobSeq;
			//				vector<vector<int>> tempJobSplit;
			//				vector <vector<vector<int>>>tempSpeedSelection;

			//				int tempMS;
			//				int tempTEC;
			//				tempJobSeq = Population[i].pJobSeq;
			//				tempJobSplit = Population[i].pJobSplit;
			//				tempMS = Population[i].MS;
			//				tempTEC = Population[i].TEC;
			//				tempSpeedSelection = Population[i].pSpeedSelection;
			//				Population[i].pJobSeq = Population[Population[i].IndexNeighbor[j]].pJobSeq;
			//				Population[i].pJobSplit = Population[Population[i].IndexNeighbor[j]].pJobSplit;
			//				Population[i].pSpeedSelection = Population[Population[i].IndexNeighbor[j]].pSpeedSelection;
			//				Population[i].MS = Population[Population[i].IndexNeighbor[j]].MS;
			//				Population[i].TEC = Population[Population[i].IndexNeighbor[j]].TEC;
			//				Population[i].getFitnessValue();
			//				Population[i].NS = 1;
			//				Population[i].pAge = 0;

			//				Population[Population[i].IndexNeighbor[j]].pJobSeq = tempJobSeq;
			//				Population[Population[i].IndexNeighbor[j]].pJobSplit = tempJobSplit;
			//				Population[Population[i].IndexNeighbor[j]].pSpeedSelection = tempSpeedSelection;
			//				Population[Population[i].IndexNeighbor[j]].MS = tempMS;
			//				Population[Population[i].IndexNeighbor[j]].TEC = tempTEC;
			//				Population[Population[i].IndexNeighbor[j]].getFitnessValue();
			//				Population[Population[i].IndexNeighbor[j]].NS = 1;

			//				Population[Population[i].IndexNeighbor[j]].pAge = 0;

			//				break;
			//				

			//			}
			//		}

			//		if (!IsChanged)
			//		{
			//			int a = rand() % PopulationSize;
			//			Population[i].pJobSeq = Population[a].pJobSeq;
			//			Population[i].pJobSplit = Population[a].pJobSplit;
			//			Population[i].pSpeedSelection = Population[a].pSpeedSelection;
			//			Population[i].MS = Population[a].MS;
			//			Population[i].TEC = Population[a].TEC;
			//			Population[i].getFitnessValue();
			//			Population[i].pAge = 0;
			//			Population[i].NS = 1;
			//		}

			//	}
			/*}*/
			if (RestartStrategy == 2)
			{

				if (Population[i].pAge > AgeTrail)
				{
					Population[i].CurrentWeightIndex++;
					bool IsChanged = false;
					for (int j = 0; j < NeighborhoodSize; j++)
					{
						if (Population[i].CurrentWeightIndex < NeighborhoodSize)
						{
							IsChanged = true;
							/*Individual temp2(1);
							Population[i] = temp2;*/
							//insert
							int pt1, pt2, pt;
							//vector<int> tempJobSeq;
							//vector<vector<int>> tempJobSplit;
							vector<vector<vector<int>>> tempspeed;

							do
							{
								pt1 = rand() % PopulationSize;
								pt2 = rand() % PopulationSize;
							} while (pt1 == pt2);
							if (pt1 < pt2)
							{
								tempspeed = Population[pt1].pSpeedSelection;
								for (pt = pt1;pt < pt2;pt++)
								{
									Population[pt].pSpeedSelection = Population[pt + 1].pSpeedSelection;
								}
								Population[pt2].pSpeedSelection = tempspeed;
							}
							else
							{
								tempspeed = Population[pt1].pSpeedSelection;
								for (pt = pt1;pt > pt2;pt--)
								{
									Population[pt].pSpeedSelection = Population[pt - 1].pSpeedSelection;
								}
								Population[pt2].pSpeedSelection = tempspeed;
							}
							break;
						}
					}

				}
			}
//			if (RestartStrategy == 3)
//			{
//				if (Population[i].pAge > AgeTrail)
//				{
//					int j = 0;
//					bool IsChanged = false;
//					while (!IsChanged)
//					{
//						if (j == i)
//						{
//							j += 1;
//						}
//
//						if (j == PopulationSize)
//						{
//							Individual temp2(1);
//							Population[i] = temp2;
//							IsChanged = true;
//							break;
//						}
//
//
//						if (Population[j].pAge > AgeTrail)
//						{
//							if (Population[i] > Population[j])
//							{
//								//交换他们的解
//								vector<int> tempJobSeq;
//								vector<vector<int>> tempJobSplit;
//								vector < vector<vector<int>>>tempSpeedSelection;
//								int tempMS;
//								int tempTEC;
//
//								tempJobSeq = Population[i].pJobSeq;
//								tempJobSplit = Population[i].pJobSplit;
//								tempMS = Population[i].MS;
//								tempTEC = Population[i].TEC;
//								tempSpeedSelection = Population[i].pSpeedSelection;
//								Population[i].pJobSeq = Population[j].pJobSeq;
//								Population[i].pJobSplit = Population[j].pJobSplit;
//								Population[i].pSpeedSelection = Population[j].pSpeedSelection;
//								Population[i].MS = Population[j].MS;
//								Population[i].TEC = Population[j].TEC;
//								Population[i].getFitnessValue();
//								Population[i].NS = 1;
//								Population[i].NSCF = 0;
//								Population[i].pAge = 0;
//
//								Population[j].pJobSeq = tempJobSeq;
//								Population[j].pJobSplit = tempJobSplit;
//								Population[j].pSpeedSelection = tempSpeedSelection;
//								Population[j].MS = tempMS;
//								Population[j].TEC = tempTEC;
//								Population[j].getFitnessValue();
//								Population[j].NS = 1;
//								Population[j].NSCF = 0;
//								Population[j].pAge = 0;
//
//								IsChanged = true;
//							}
//						}
//
//						j++;
//
//						/*if (j == PopulationSize - 1)
//						{
//						Individual temp2(2);
//						Population[i] = temp2;
//						IsChanged = true;
//						}*/
//					}
//
//				}
//			}
//			if (RestartStrategy == 4)
//			{
//				vector<int> tempIndex;
//				for (int i = 0; i < PopulationSize; i++)
//				{
//					if (Population[i].pAge > AgeTrail)
//					{
//						tempIndex.push_back(i);
//					}
//				}
//
//				if (tempIndex.size() == 1)
//				{
//					Individual temp2(1);
//					Population[tempIndex[0]] = temp2;
//					break;
//				}
//				if (tempIndex.size() > 1)
//				{
//					for (int i = 0; i < tempIndex.size() - 1; i++)
//					{
//						//完成两个解的互换
//						//交换他们的解
//						vector<int> tempJobSeq;
//						vector<vector<int>> tempJobSplit;
//						vector < vector<vector<int>>>tempSpeedSelection;
//						int tempMS;
//						int tempTEC;
//
//						tempJobSeq = Population[tempIndex[i]].pJobSeq;
//						tempJobSplit = Population[tempIndex[i]].pJobSplit;
//						tempSpeedSelection = Population[i].pSpeedSelection;
//						tempMS = Population[tempIndex[i]].MS;
//						tempTEC = Population[tempIndex[i]].TEC;
//
//						Population[tempIndex[i]].pJobSeq = Population[tempIndex[i + 1]].pJobSeq;
//						Population[tempIndex[i]].pJobSplit = Population[tempIndex[i + 1]].pJobSplit;
//						Population[tempIndex[i]].pSpeedSelection = Population[tempIndex[i + 1]].pSpeedSelection;
//						Population[tempIndex[i]].MS = Population[tempIndex[i + 1]].MS;
//						Population[tempIndex[i]].TEC = Population[tempIndex[i + 1]].TEC;
//						Population[tempIndex[i]].getFitnessValue();
//						Population[tempIndex[i]].NS = 1;
//						Population[tempIndex[i]].NSCF = 0;
//						Population[tempIndex[i]].pAge = 0;
//
//						Population[i + 1].pJobSeq = tempJobSeq;
//						Population[i + 1].pJobSplit = tempJobSplit;
//						Population[i + 1].pSpeedSelection = tempSpeedSelection;
//						Population[i + 1].MS = tempMS;
//						Population[i + 1].TEC = tempTEC;
//						Population[i + 1].getFitnessValue();
//						Population[i + 1].NS = 1;
//						Population[i + 1].NSCF = 0;
//						Population[i + 1].pAge = 0;
//					}
//				}
//			}
//

		}
	}
}

void MDABC::NondominatedSorting()
{
	//re-initialize the n value
	for (int i = 0; i < PopulationSize; i++)
	{
		Population[i].n = 0;
		Population[i].distance = 0;
	}
	int maxN = 0;
	for (int i = 0; i < PopulationSize; i++)
	{
		for (int j = 0; j < PopulationSize; j++)
		{
			if (Population[i] < Population[j])
			{
				Population[i].n++;
			}
		}
		if (Population[i].n > maxN)
		{
			maxN = Population[i].n;
		}
	}
	//////generate nondominated fronts method(1)
	vector<vector<Individual>> p;
	p.resize(maxN + 1);
	for (int i = 0; i < PopulationSize; i++)
	{
		p[Population[i].n].push_back(Population[i]);
	}
	////calculate the crowding distance in P[i]
	//for (int i = 0; i < p.size(); i++)
	//{
	//	if (p[i].size() > 0)
	//	{
	//		//compare MS
	//		sort(p[i].begin(), p[i].end(), CompareMS());
	//		//给边界点之外的点赋值
	//		for (int j = 1; j < p[i].size() - 1; j++)
	//		{
	//			p[i][j].distance += abs(p[i][j + 1].MS - p[i][j - 1].MS);
	//		}
	//		//给边界点赋值 赋值大值保证能被选中
	//		p[i][0].distance += 10000;
	//		p[i][p[i].size() - 1].distance += 10000;
	//		//compare TEC
	//		sort(p[i].begin(), p[i].end(), CompareTEC());
	//		//给边界点之外的点赋值
	//		for (int j = 1; j < p[i].size() - 1; j++)
	//		{
	//			p[i][j].distance += abs(p[i][j + 1].TEC - p[i][j - 1].TEC);
	//		}
	//		//给边界点赋值 赋值大值保证能被选中
	//		p[i][0].distance += 10000;
	//		p[i][p[i].size() - 1].distance += 10000;


	//		sort(p[i].begin(), p[i].end(), CompareDistance());
	//	}
	//}
	////update the population
	//Population.clear();
	//bool isFull = false;
	//for (int i = 0; i < p.size(); i++)
	//{
	//	if (!isFull)
	//	{
	//		for (int j = 0; j < p[i].size(); j++)
	//		{
	//			Population.push_back(p[i][j]);
	//			if (Population.size() >= PopulationSize)
	//			{
	//				isFull = true;
	//				break;
	//			}
	//		}
	//	}
	//	
	//}
	////update ParetoPopulation
	//UpdateParetoPouplation();
}

void MDABC::UpdateParetoPouplation(Individual& Ind)
{
	/*if (Ind.MS < RP.MS)
	{
		RP.MS = Ind.MS;
	}
	if (Ind.TEC < RP.TEC)
	{
		RP.TEC = Ind.TEC;
	}*/
	////remove from ParetoPopulation dominated by this
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

void MDABC::OutputParetoPopulation(int ins, double costTime)
{
	sort(ParetoPopulation.begin(), ParetoPopulation.end(), CompareMS());
	/*for (int i = 0; i < ParetoPopulation.size(); i++)
	{
	cout << ParetoPopulation[i].MS << " " << ParetoPopulation[i].MF << " " << ParetoPopulation[i].TF << endl;
	}*/

	char outfile[40];
	sprintf_s(outfile, "results.mdabc\\%d_%d_%d_%d.txt", pJob, pStage, pType, ins);
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
