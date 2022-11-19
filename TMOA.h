#pragma once
#include "Individual.h"
#include "ParetoPoint.h"
#include "CompareIndividual.h"
#include <iomanip>
class TMOA
{
public:
	TMOA();
	TMOA(int _H, int _NeighborhoodSize, int _NeighborhoodTrial, int _AgeTrial, int _MaximumCPU);
	~TMOA();
public:
	void Gen2WV_NBI(); //generate weight vectors for two objectives
	void UpdateParetoPouplation(Individual& Ind);
	void UpdateParetoPouplation();
	void InitializePopulation();  //initialize the population
	void EvolvePopulation(long time);      //Evolve the population
	void NondominatedSorting();     //non-dominated sorting
	//void ConductImprovementStrategy();

	void run(int ins);


public:
	vector<Individual>  Population;   //current population
	vector<Individual>  ParetoPopulation;   //pareto points
public:
	int PopulationSize;    //the current population size
	int NeighborhoodSize;  //the neighborhood size
	int MaximumCPUTime;  //the maximum CPU time
	int AgeTrail;          //age trial
	int NeighborhoodTrial;   //neighborhood trial

public:
	//for the weight vectors of two objectives
	int H;
	vector<double> W1;
	vector<double> W2;

public:
	void TwoPointsXoverAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd);
	void UpdateReferencePoint(Individual& Ind);
	void OutputParetoPopulation(int ins, double costTime);

public:
	ParetoPoint RP;  //reference point 
};
TMOA::TMOA()
{

}

TMOA::~TMOA()
{

}
TMOA::TMOA(int _H, int _NeighborhoodSize, int _NeighborhoodTrial, int _AgeTrial, int _MaximumCPU)
{
	this->H = _H;
	this->NeighborhoodSize = _NeighborhoodSize;
	this->NeighborhoodTrial = _NeighborhoodTrial;
	this->AgeTrail = _AgeTrial;
	this->MaximumCPUTime = _MaximumCPU * pJob * pStage;

	this->PopulationSize = H + 1; //149->150  199->200  249->250
}
void TMOA::Gen2WV_NBI()
{
	W1.resize(PopulationSize);
	W2.resize(PopulationSize);


	for (int i = 0, j = H; i <= H, j >= 0; i++, j--)
	{
		W1[i] = i * 1.0 / H;
		W2[i] = j * 1.0 / H;
	}
}
void TMOA::run(int ins)
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
void TMOA::InitializePopulation()
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
	//get the fitness value
	for (int i = 0; i < PopulationSize; i++)
	{
		Population[i].getFitnessValue();
	}
}
void TMOA::UpdateReferencePoint(Individual& Ind)
{
	if (Ind.MS < RP.MS)
	{
		RP.MS = Ind.MS;
	}
	if (Ind.TEC < RP.TEC)
	{
		RP.TEC = Ind.TEC;
	}
}
void TMOA::EvolvePopulation( long time)
{
	//int generations = 0;
	long InitTime = GetTickCount();
	while ((GetTickCount() - InitTime) < MaximumCPUTime)
		//while(generations < 1000)
	{
		//generations++;
		//the first phase
		for (int i = 0; i < PopulationSize; i++)
		{
			Individual Neighbor;
			Neighbor = Population[i];
			Population[i].getNeigbors(Neighbor, Population[i].CurrentNeighbrohoodIndex);
			//Individual* Neighbor = Population[i].getNeigbors(Population[i].CurrentNeighbrohoodIndex);
			UpdateParetoPouplation(Neighbor);
			//update the current population
			if (Neighbor.Fitness < Population[i].Fitness)//将适应度越小的
			{
				Population[i].pJobSeq = Neighbor.pJobSeq;
				Population[i].pJobSplit = Neighbor.pJobSplit;
				Population[i].pSpeedSelection = Neighbor.pSpeedSelection;
				Population[i].MS = Neighbor.MS;
				Population[i].TEC = Neighbor.TEC;
				Population[i].CurrentNeighbrohoodIndex = 1;
				Population[i].CurrentWeightIndex = 0;
				Population[i].getFitnessValue();
				Population[i].pAge = 0;
				Population[i].neighborhoodAge = 0;
			}
			else
			{
				Population[i].pAge++;
				Population[i].neighborhoodAge++;

				if (Population[i].neighborhoodAge > NeighborhoodTrial)
				{
					Population[i].CurrentNeighbrohoodIndex++;
				}

				if (Population[i].CurrentNeighbrohoodIndex > 8)
				{

					Population[i].CurrentNeighbrohoodIndex = 1;
				}
			}

		}
		//NondominatedSorting();
		// Collaborative global search phase
		//evolution
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


			Individual tempInd;

			int r3 = rand() % NeighborhoodSize;
			int a = Population[maxIndex].IndexNeighbor[r3];

			TwoPointsXoverAndMutation(Population[maxIndex], Population[a], tempInd);

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
						Population[select[k]].pSpeedSelection = tempInd.pSpeedSelection;
						Population[select[k]].pJobSplit = tempInd.pJobSplit;
						Population[select[k]].MS = tempInd.MS;
						Population[select[k]].TEC = tempInd.TEC;
						Population[select[k]].Fitness = tempInd.Fitness;
						Population[select[k]].CurrentNeighbrohoodIndex = 1;
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
							Population[select[k]].pJobSplit = tempInd.pJobSplit;
							Population[select[k]].pSpeedSelection = tempInd.pSpeedSelection;
							Population[select[k]].MS = tempInd.MS;
							Population[select[k]].TEC = tempInd.TEC;
							Population[select[k]].CurrentWeightIndex = 0;
							Population[select[k]].getFitnessValue();
							Population[select[k]].CurrentNeighbrohoodIndex = 1;
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

		//the third phase
		for (int i = 0; i < PopulationSize; i++)
		{
			if (Population[i].pAge > AgeTrail)
			{
				Population[i].CurrentWeightIndex++;

				if (Population[i].CurrentWeightIndex < NeighborhoodSize)
				{
					Population[i].getFitnessValue();
					Population[i].pAge = 0;
				}
				else
				{
					//cout << Population[i].CurrentWeightIndex << "A loop is finished!" << endl;
					//cout << "A loop is finished!" << endl;
					Population[i].CurrentWeightIndex = 0;
					Population[i].pAge = 0;

					//conduct blocking swapping operator
					//Population[i].conductPermutation();

					Population[i].getObjectives2(1);
					Population[i].getFitnessValue();
				}

			}
		}
	}
}
void  TMOA::TwoPointsXoverAndMutation(Individual& Ind1, Individual& Ind2, Individual& ResultInd)
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
void  TMOA::UpdateParetoPouplation(Individual& Ind)
{
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
//void  TMOA::UpdateParetoPopulation()
//{
//	//conduct the energy saving procedure
//	/*for (int i = 0; i < ParetoPopulation.size(); i++)
//	{
//		ParetoPopulation[i].getObjectives4(2);
//	}*/
//	vector<Individual> tempPopulation = ParetoPopulation;
//
//	ParetoPopulation.clear();
//
//	//update the pareto population
//	for (int i = 0; i < tempPopulation.size(); i++)
//	{
//		//remove from ParetoPopulation dominated by this
//		for (int n = 0; n < ParetoPopulation.size(); n++)
//		{
//			if (tempPopulation[i] > ParetoPopulation[n])
//			{
//				ParetoPopulation.erase(ParetoPopulation.begin() + n);
//				n--;
//			}
//		}
//		//Add this to ParetoPopulation if no points dominate this
//		bool  isDominated = false;
//		for (vector<Individual>::iterator it = ParetoPopulation.begin(); it != ParetoPopulation.end(); it++)
//		{
//			if (*it > tempPopulation[i] || *it == tempPopulation[i])
//			{
//				isDominated = true;
//				break;
//			}
//		}
//
//		if (!isDominated)
//		{
//			ParetoPopulation.push_back(tempPopulation[i]);
//		}
//	}
//}
void TMOA::OutputParetoPopulation(int ins, double costTime)
{
	//update the paretopopulation



	sort(ParetoPopulation.begin(), ParetoPopulation.end(), CompareMS());
	/*for (int i = 0; i < ParetoPopulation.size(); i++)
	{
	cout << ParetoPopulation[i].MS << " " << ParetoPopulation[i].MF << " " << ParetoPopulation[i].TF << endl;
	}*/

	char outfile[40];
	sprintf_s(outfile, "results.tmoa\\%d_%d_%d_%d.txt", pJob, pStage, pType, ins);
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