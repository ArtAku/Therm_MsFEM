#pragma once
#include "libraries.h"

struct boundaryData1D
{
	std::vector<double> _nodeValue;
	std::vector<int> _finElemNode[2];
	double _boundaryValue[2];
	std::vector<double> _lambda;
	std::vector<int> _nodeIndex;
};

boundaryData1D insertNode(boundaryData1D data, double nodeValue, int nodeIndex)
{
	boundaryData1D result = data;

	std::vector<double>::iterator it = result._nodeValue.begin();
	std::vector<int>::iterator it1 = result._nodeIndex.begin();
	while (it != result._nodeValue.end())
	{
		if (nodeValue < *it)
		{
			result._nodeValue.insert(it, nodeValue);
			result._nodeIndex.insert(it1, nodeIndex);
			break;
		}
		it++;
		it1++;
	}
	if (data._nodeValue.size() == result._nodeValue.size())
	{
		result._nodeValue.push_back(nodeValue);
		result._nodeIndex.push_back(nodeIndex);
	}
	return result;
}