#include "two_dimensional_task.h"

//struct boundaryData1D
//{
//	std::vector<double> _nodeValue;
//	std::vector<int> _finElemNode[2];
//	double _boundaryValue[2];
//	std::vector<double> _lambda;
//	std::vector<int> _nodeIndex;
//};
//
//boundaryData1D insertNode(boundaryData1D data, double nodeValue, int nodeIndex)
//{
//	boundaryData1D result = data;
//
//	std::vector<double>::iterator it = result._nodeValue.begin();
//	std::vector<int>::iterator it1 = result._nodeIndex.begin();
//	while (it != result._nodeValue.end())
//	{
//		if (nodeValue < *it)
//		{
//			result._nodeValue.insert(it, nodeValue);
//			result._nodeIndex.insert(it1, nodeIndex);
//			break;
//		}
//		it++;
//		it1++;
//	}
//	if (data._nodeValue.size() == result._nodeValue.size())
//	{
//		result._nodeValue.push_back(nodeValue);
//		result._nodeIndex.push_back(nodeIndex);
//	}
//	return result;
//}

two_dimensional_task::two_dimensional_task(std::map<int, int> nodeIndex, std::vector<double> pointX, std::vector<double> pointY, std::vector<int> finElemNode[3], std::vector<double> lambda, std::vector<int> boundaryNode, int basisNode)
{
	_G = new double*[3];
	for (int i(0); i < 3; i++)
		_G[i] = new double[3];
	_finElemNode = new std::vector<int>[3];

	_n = pointX.size();
	_nodeIndex = nodeIndex;
	_pointX = pointX;
	_pointY = pointY;
	for (int i(0); i < 3; i++)
		_finElemNode[i] = finElemNode[i];
	_lambda = lambda;
	_finElemCount = finElemNode[0].size();
	_boundaryNode = boundaryNode;
	_basisNode = basisNode;

	_b.reserve(_n);
	_x.reserve(_n);
	_ig.reserve(_n + 1);

	genPortrait();
}

void two_dimensional_task::compileSLAU()
{
	int curr_node, v, count;
	for (int i(0); i < _signMatrixElemCount; i++)
		_gg.push_back(0);
	for (int i(0); i < _n; i++)
		_b.push_back(0);
	for (int i(0); i < (int)_finElemNode[0].size(); i++)
	{
		calcLocalMatrix(i);
		for (int j(0); j < 3; j++)
		{
			curr_node = _finElemNode[j][i];
			for (int k(0); k < 3; k++)
			{
				count = 0;
				for (auto var : _sortNode[curr_node])
				{
					if (var == _finElemNode[k][i])
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

void two_dimensional_task::imposeBoundaryCond()
{
	int place;
	int count;
	std::set<int>::iterator it;
	std::map<int, int> nodeIndex[4];
	boundaryData1D lineData[4];
	double one_dim_boundary_nodes[4]{ 0 , 0, 0, 0 };
	one_dim_boundary_nodes[_basisNode] = 1;
	std::vector<double>::iterator minX = std::min_element(_pointX.begin(), _pointX.end());
	std::vector<double>::iterator maxX = std::max_element(_pointX.begin(), _pointX.end());
	std::vector<double>::iterator minY = std::min_element(_pointY.begin(), _pointY.end());
	std::vector<double>::iterator maxY = std::max_element(_pointY.begin(), _pointY.end());
	for (auto bound : _boundaryNode)
	{
		if (_pointX[_nodeIndex[bound]] == *minX)
			lineData[0] = insertNode(lineData[0], _pointY[_nodeIndex[bound]], bound);
		if (_pointY[_nodeIndex[bound]] == *maxY)
			lineData[1] = insertNode(lineData[1], _pointX[_nodeIndex[bound]], bound);
		if (_pointX[_nodeIndex[bound]] == *maxX)
			lineData[2] = insertNode(lineData[2], _pointY[_nodeIndex[bound]], bound);
		if (_pointY[_nodeIndex[bound]] == *minY)
			lineData[3] = insertNode(lineData[3], _pointX[_nodeIndex[bound]], bound);
	}
	for (int i(0); i < 4; i++)
	{
		for (int j(0); j < (int)lineData[i]._nodeIndex.size(); j++)
			nodeIndex[i].insert(std::pair<int, int>(lineData[i]._nodeIndex[j], j));
		lineData[i]._boundaryValue[0] = one_dim_boundary_nodes[i];
		place = i < 3 ? i + 1 : 0;
		lineData[i]._boundaryValue[1] = one_dim_boundary_nodes[place];

		for (int j(0); j < (int)lineData[i]._nodeIndex.size() - 1; j++)
		{
			lineData[i]._finElemNode[0].push_back(lineData[i]._nodeIndex[j]);
			lineData[i]._finElemNode[1].push_back(lineData[i]._nodeIndex[j + 1]);
			for (int k(0); k < (int)_finElemNode[0].size(); k++)
			{
				if ((lineData[i]._finElemNode[0][j] == _finElemNode[0][k] || lineData[i]._finElemNode[0][j] == _finElemNode[1][k] || lineData[i]._finElemNode[0][j] == _finElemNode[2][k])
					&& (lineData[i]._finElemNode[1][j] == _finElemNode[0][k] || lineData[i]._finElemNode[1][j] == _finElemNode[1][k] || lineData[i]._finElemNode[1][j] == _finElemNode[2][k]))
				{
					lineData[i]._lambda.push_back(_lambda[k]);
					break;
				}
			}
		}
	}

	std::vector<double> res[4];

	one_dimensional_task left_bound(nodeIndex[0], lineData[0]._nodeValue, lineData[0]._finElemNode, lineData[0]._lambda, lineData[0]._boundaryValue);
	one_dimensional_task top_bound(nodeIndex[1], lineData[1]._nodeValue, lineData[1]._finElemNode, lineData[1]._lambda, lineData[1]._boundaryValue);
	one_dimensional_task right_bound(nodeIndex[2], lineData[2]._nodeValue, lineData[2]._finElemNode, lineData[2]._lambda, lineData[2]._boundaryValue);
	one_dimensional_task bottom_bound(nodeIndex[3], lineData[3]._nodeValue, lineData[3]._finElemNode, lineData[3]._lambda, lineData[3]._boundaryValue);
	res[0] = left_bound.runSolver();
	res[1] = top_bound.runSolver();
	res[2] = right_bound.runSolver();
	res[3] = bottom_bound.runSolver();

	for (int i(0); i < 4; i++)
	{
		for (int k(0); k < (int)lineData[i]._nodeIndex.size(); k++)
		{
			count = 0;
			it = _sortNode[lineData[i]._nodeIndex[k]].find(lineData[i]._nodeIndex[k]);
			for (std::set<int>::iterator j = _sortNode[lineData[i]._nodeIndex[k]].begin(); j != _sortNode[lineData[i]._nodeIndex[k]].end(); j++, count++)
				if (it == j)
					_gg[_ig[lineData[i]._nodeIndex[k]] + count] = 1;
				else
					_gg[_ig[lineData[i]._nodeIndex[k]] + count] = 0;
			_b[lineData[i]._nodeIndex[k]] = res[i][k];
		}
	}
}

void two_dimensional_task::calcLocalMatrix(int cur_el)
{
	double detD = (_pointX[_finElemNode[1][cur_el]] - _pointX[_finElemNode[0][cur_el]]) * (_pointY[_finElemNode[2][cur_el]] - _pointY[_finElemNode[0][cur_el]]) - (_pointX[_finElemNode[2][cur_el]] - _pointX[_finElemNode[0][cur_el]]) * (_pointY[_finElemNode[1][cur_el]] - _pointY[_finElemNode[0][cur_el]]);

	_alpha[0][0] = (_pointX[_finElemNode[1][cur_el]] * _pointY[_finElemNode[2][cur_el]] - _pointX[_finElemNode[2][cur_el]] * _pointY[_finElemNode[1][cur_el]]) / detD;
	_alpha[0][1] = (_pointY[_finElemNode[1][cur_el]] - _pointY[_finElemNode[2][cur_el]]) / detD;
	_alpha[0][2] = (_pointX[_finElemNode[2][cur_el]] - _pointX[_finElemNode[1][cur_el]]) / detD;

	_alpha[1][0] = (_pointX[_finElemNode[2][cur_el]] * _pointY[_finElemNode[0][cur_el]] - _pointX[_finElemNode[0][cur_el]] * _pointY[_finElemNode[2][cur_el]]) / detD;
	_alpha[1][1] = (_pointY[_finElemNode[2][cur_el]] - _pointY[_finElemNode[0][cur_el]]) / detD;
	_alpha[1][2] = (_pointX[_finElemNode[0][cur_el]] - _pointX[_finElemNode[2][cur_el]]) / detD;

	_alpha[2][0] = (_pointX[_finElemNode[0][cur_el]] * _pointY[_finElemNode[1][cur_el]] - _pointX[_finElemNode[1][cur_el]] * _pointY[_finElemNode[0][cur_el]]) / detD;
	_alpha[2][1] = (_pointY[_finElemNode[0][cur_el]] - _pointY[_finElemNode[1][cur_el]]) / detD;
	_alpha[2][2] = (_pointX[_finElemNode[1][cur_el]] - _pointX[_finElemNode[0][cur_el]]) / detD;

	for (int i(0); i < 3; i++)
		for (int j(0); j < 3; j++)
			_G[i][j] = _lambda[cur_el] * abs(detD) * (_alpha[i][1] * _alpha[j][1] + _alpha[i][2] * _alpha[j][2]) / 2;
}

std::vector<double> two_dimensional_task::runSolver()
{
	compileSLAU();
	imposeBoundaryCond();
	local_optimal_scheme solverSLAU(_n, _signMatrixElemCount, 1e-20, _gg, _b, _ig, _jg);
	_x = solverSLAU.getSolve();
	return _x;

}

void two_dimensional_task::writeFile()
{
	std::fstream output("output.xls");
	double t_cur = 50.0;
	double t_step = 3;
	double curPoint[2] = { 0.0 , 1.0 };
	while (curPoint[1] >= 0)
	{
		while (curPoint[0] <= 1)
		{
			output << getSolutionInPoint(curPoint[0], curPoint[1]) << "\t";
			curPoint[0] += 0.1;
		}
		output << std::endl;
		curPoint[1] -= 0.1;
		curPoint[0] = 0;
		t_cur -= t_step;
	}
}

double two_dimensional_task::getSolutionInPoint(double x, double y)
{
	double ans = 0;
	double check[3];
	int elemInd;

	for (int i(0); i < (int)_finElemNode[0].size(); i++)
	{
		check[0] = (_pointX[_finElemNode[0][i]] - x) * (_pointY[_finElemNode[1][i]] - _pointY[_finElemNode[0][i]]) - (_pointX[_finElemNode[1][i]] - _pointX[_finElemNode[0][i]]) * (_pointY[_finElemNode[0][i]] - y);
		check[1] = (_pointX[_finElemNode[1][i]] - x) * (_pointY[_finElemNode[2][i]] - _pointY[_finElemNode[1][i]]) - (_pointX[_finElemNode[2][i]] - _pointX[_finElemNode[1][i]]) * (_pointY[_finElemNode[1][i]] - y);
		check[2] = (_pointX[_finElemNode[2][i]] - x) * (_pointY[_finElemNode[0][i]] - _pointY[_finElemNode[2][i]]) - (_pointX[_finElemNode[0][i]] - _pointX[_finElemNode[2][i]]) * (_pointY[_finElemNode[2][i]] - y);
		if (!(check[0] * check[1] < 0 || check[1] * check[2] < 0 || check[0] * check[2] < 0))
		{
			elemInd = i;
			break;
		}
	}
	calcLocalMatrix(elemInd);
	for (int i(0); i < 3; i++)
		ans += (_alpha[i][0] + _alpha[i][1] * x + _alpha[i][2] * y) * _x[_finElemNode[i][elemInd]];
	return ans;
}
