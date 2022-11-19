#pragma once
#include "Individual.h"
class CompareMS
{
public:
	bool operator()(Individual& obj1, Individual& obj2)
	{
		return obj1.MS < obj2.MS;
	}
};


class CompareTEC
{
public:
	bool operator()(Individual& obj1, Individual& obj2)
	{
		return obj1.TEC < obj2.TEC;
	}
};
class CompareDistance
{
public:
	bool operator()(Individual& obj1, Individual& obj2)
	{
		return obj1.distance > obj2.distance;
	}
};