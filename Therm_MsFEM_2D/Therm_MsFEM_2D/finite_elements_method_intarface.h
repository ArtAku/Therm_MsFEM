#pragma once
#include "libraries.h"
class finite_elements_method_intarface
{
public:	
	virtual std::vector<double> runSolver() = 0; // Get solve of a task. Returns _x;

	int _signMatrixElemCount; // Number of matrix significant elements.
	int _n; // Count of nodes.
	int _finElemCount; // Count of finite elements.
	
protected:
	virtual void genPortrait()
	{
		int finElemCount = _finElemNode[0].size();
		_sortNode = new std::set<int>[_n];
		for (int i(0); i < finElemCount; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					_sortNode[_nodeIndex.at(_finElemNode[j][i])].insert(_nodeIndex.at(_finElemNode[k][i]));

		int signMatrixElemCount = 0;
		for (int i(0); i < _n; i++)
			_signMatrixElemCount = _signMatrixElemCount + _sortNode[i].size();
		_gg.reserve(_signMatrixElemCount);
		_jg.reserve(_signMatrixElemCount);

		_ig.push_back(0);
		for (int i(1), j(0); i <= _n; i++)
		{
			j += _sortNode[i - 1].size();
			_ig.push_back(j);
		}

		for (int i(0), j(0); i < _n; i++)
			for (auto var : _sortNode[i])
				_jg.push_back(var);
	}; // Fill ig and jg (generate matrix portrait).
	virtual void calcLocalMatrix(int cur_el) = 0; // Caclulate local matrix G for finite element cur_el.
	virtual void compileSLAU() = 0; // Create a global matrix of local matrices.
	virtual void imposeBoundaryCond() = 0; // Impose boundary conditions on SLAU.

	// Matrix vectors 
	std::vector<double> _gg; // Values of each matrix element. Size = _signMatrixElemCount.
	std::vector<int> _ig; // Index of first element of each line. Size = _n + 1.
	std::vector<int> _jg; // Column number of each matrix element. Size = _signMatrixElemCount.

	std::vector<double> _x; // desired vector of SLAU. Size = _n.
	std::vector<double> _b; // Right-part vector of SLAU. Size = _n.

	double** _G; // local matrix for finite element. Size related with basis functions and dimension of task.

	std::vector<double> _lambda; // thermal coefficient. Size = _finElemCount.
	std::vector<int>* _finElemNode; // Nodes for finiteElements. Size1 = related with basis functions and dimension of task. Size2 = _finElemCount.

	std::set<int>* _sortNode; // list of neighbour-nodes. Size1 = _n. Size2 = related with basis functions and dimension of task.
	std::map<int, int> _nodeIndex; // Pair of relations global node - local node. Size = _n;

};

