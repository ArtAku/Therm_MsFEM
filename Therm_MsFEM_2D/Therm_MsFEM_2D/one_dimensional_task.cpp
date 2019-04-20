#include "one_dimensional_task.h"



one_dimensional_task::one_dimensional_task(std::map<int, int> nodeIndex, std::vector<double> point, std::vector<int> finElemNode[2], std::vector<double> lambda, double boundaryX[2])
{
	_G = new double*[2];
	for (int i(0); i < 2; i++)
		_G[i] = new double[2];
	_finElemNode = new std::vector<int>[2];

	_point = point;
	_nodeIndex = nodeIndex;
	for (int i(0); i < 2; i++)
	{
		_finElemNode[i] = finElemNode[i];
		_boundaryX[i] = boundaryX[i];
	}
	_lambda = lambda;
	_n = point.size();
	_finElemCount = finElemNode[0].size();
	// Allocating memory
	_ig.reserve(_n + 1);
	_b.reserve(_n);
	_x.reserve(_n);

	
}

void one_dimensional_task::imposeBoundaryCond()
{
	// First node
	_gg[0] = 1;
	_gg[1] = 0;
	_b[0] = _boundaryX[0];

	// Last node
	_gg[_gg.size() - 1] = 1;
	_gg[_gg.size() - 2] = 0;
	_b[_b.size() - 1] = _boundaryX[1];
}

void one_dimensional_task::calcLocalMatrix(int cur_el) 
{
	_G[0][0] = _lambda[cur_el] / abs(_point[_nodeIndex.at(_finElemNode[1][cur_el])] - _point[_nodeIndex.at(_finElemNode[0][cur_el])]) / 2;
	_G[0][1] = -_G[0][0];
	_G[1][0] = -_G[0][0];
	_G[1][1] = _G[0][0];
}

void one_dimensional_task::compileSLAU()
{
	int curr_node, v, count;
	// Initialize arrays of values with 0
	for (int i(0); i < _signMatrixElemCount; i++)
		_gg.push_back(0);
	for (int i(0); i < _n; i++)
		_b.push_back(0);

	// Impose every finite element on global matrix
	for (int i(0); i < _finElemCount; i++)
	{
		calcLocalMatrix(i);
		for (int j(0); j < 2; j++)
		{
			curr_node = _nodeIndex.at(_finElemNode[j][i]);
			for (int k(0); k < 2; k++)
			{
				count = 0;
				for (auto var : _sortNode[curr_node])
				{
					if (var == _nodeIndex.at(_finElemNode[k][i]))
					{
						v = count;
						break;
					}
					count++;
				}
				_gg[_ig[curr_node] + v] += _G[j][k];
			}
		}
	}
}

std::vector<double> one_dimensional_task::runSolver()
{
	genPortrait();
	compileSLAU();
	imposeBoundaryCond();
	local_optimal_scheme solverSLAU(_n, _signMatrixElemCount, 1e-20, _gg, _b, _ig, _jg);
	_x = solverSLAU.getSolve();
	return _x;
}
