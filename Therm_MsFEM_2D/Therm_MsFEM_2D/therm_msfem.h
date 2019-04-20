#pragma once
#include "two_dimensional_task.h"
#include "libraries.h"

class therm_msfem
{
public:
	therm_msfem(std::string in);
	
	void runSolve();

private:
	std::vector<int> _finElemNode[3];
	std::map<int, int> _nodeIndex;
	std::vector<double> _pointX;
	std::vector<double> _pointY;
	std::vector<double> _lambda;
	std::vector<int> _boundaryNode;
	int _basisNode;
};

