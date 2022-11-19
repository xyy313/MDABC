#pragma once


#include <cstdlib>
#include <cmath>
#include <fstream>
#include <Windows.h>
#include <iostream>
#include <vector>
#include <map>
#include <ctime>
#include <string>
#include <time.h>
#include <algorithm>
#include "ParetoPoint.h"

template<typename T>
struct Pair1
{
	int dim;
	T value;
	T weightMS;
	T weightTEC;
};


template<typename T>
class PairGreater1
{
public:
	bool operator () (Pair1<T> a, Pair1<T> b)
	{
		return a.value > b.value;
	}
};

template<typename T>
class PairLess1
{
public:
	bool operator () (Pair1<T> a, Pair1<T> b)
	{
		return a.value < b.value;
	}
};

