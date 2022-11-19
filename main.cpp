#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <tchar.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <numeric>
#include <Windows.h>
#include <vector>
#include "MDABC.h"
#include "MOEAD.h"
#include "NSGAII.h"
#include "TMOA.h"
#include "MOCGWO.h"
using namespace std;

void main()
{
	//�����㷨����
	//��������ʱ��
	int CPU = 100;

	int NumberofEvluations = 100;
	//������Ⱥ��ģ
	int HValue = 99;
	int PopulationSize = 99;

	//���������С
	int NeighborhoodSize = 20;

	//����VNS
	int VNSTrial = 1;

	//����ABC
	int AgeTrial = 50;
	int NeighborhoodTrial = 8;
	//���ڶ�̬�������
	bool IsVD = 1;

	//������������
	int ReStart = 1;
	int NumberOfInstances = 0;
	//����������
	int MaximumFuncEvals =200;
	for (int i = 0; i < 5; i++)  //for the number of jobs  5
	{
		for (int j = 0; j < 3; j++)  //for the number of stages  4
		{
			for (int k = 0; k < 4; k++)
			{
				for (int ins = 0; ins < 5; ins++)  //for the instances of each combination  10
				{
					NumberOfInstances++;
					int InJob = Jobs[i];
					int InStage = Stages[j];
					int InType = Types[k];
					//int InSetup = SetupTime[k];
					//GenerateInstances(InJob, InStage, InType,ins);  //output the instances
					GenerateInstances(InJob, InStage, InType, NumberOfInstances+1234);
					//GenerateInstances_my();
					GetUpperAndLowerBounds(InType);
					for (int rep = 0; rep < 5; rep++)
					{
						MDABC mdabc(HValue, NeighborhoodSize, NeighborhoodTrial,  AgeTrial, IsVD, CPU);
						mdabc.run(ins);
						//MOEAD moead(HValue, NeighborhoodSize, NeighborhoodTrial, AgeTrial, CPU);
						//moead.run(ins);
						//NSGAII nsga2(PopulationSize, NeighborhoodSize, MaximumFuncEvals);
						//nsga2.run(ins);
						//TMOA tmoa(HValue, NeighborhoodSize, NeighborhoodTrial, AgeTrial, CPU);
						//tmoa.run(ins);
						//MOCGWO mocgwo(PopulationSize, CPU);
						//mocgwo.run(ins);
						cout << pJob << "_" << pStage << "_" << pType << "_" << ins << " rep:" << rep << "End" << endl;
					}
				}
			}
		}
	}

}