#pragma once
#include "libraries.h"
#include "one_dimensional_task.h"
#include "finite_elements_method_intarface.h"
#include "boundaryData.h"

class two_dimensional_task:finite_elements_method_intarface
{
public:
	two_dimensional_task(std::map<int, int> nodeIndex, std::vector<double> pointX, std::vector<double> pointY, std::vector<int> finElemNode[3], std::vector<double> lambda, std::vector<int> boundaryNode, int basisNode);
	std::vector<double> runSolver();
	double getSolutionInPoint(double x, double y);
	void writeFile();

private:
	void calcLocalMatrix(int cur_el) override; // Caclulate local matrix G for finite element cur_el
	void compileSLAU() override; // Create a global matrix of local matrices
	void imposeBoundaryCond() override; // Impose boundary conditions on SLAU 

	std::vector<double> _pointX; // Coordinates of nodes axis x. Size = _n.
	std::vector<double> _pointY; // Coordinates of nodes axis y. Size = _n.
	double _alpha[3][3]; // Matrix for L-coordinates
	int _basisNode; // One of 4 nodes which we will initialize with 1, others with 0
	std::vector<int> _boundaryNode; // List of boundary nodes. Size = idc.
};