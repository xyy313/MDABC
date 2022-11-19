#pragma once
class ParetoPoint
{
public:
	ParetoPoint();
	ParetoPoint(int _MS, int _TEC);
	~ParetoPoint();
public:
	bool operator==(ParetoPoint& obj);		//equal
	bool operator>(ParetoPoint& obj);		//dominate
	bool operator<(ParetoPoint& obj);		//dominated
public:
	double MS;		//makeSpan
	double TEC;		//total energy cost
};


ParetoPoint::ParetoPoint()
{
}

ParetoPoint::ParetoPoint(int _MS, int _TEC)
{
	this->MS = _MS;
	this->TEC = _TEC;
}

//equal
bool ParetoPoint::operator==(ParetoPoint& obj)
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

//dominate
bool ParetoPoint::operator>(ParetoPoint& obj)
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

//dominated
bool ParetoPoint::operator<(ParetoPoint& obj)
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
ParetoPoint::~ParetoPoint()
{
}


