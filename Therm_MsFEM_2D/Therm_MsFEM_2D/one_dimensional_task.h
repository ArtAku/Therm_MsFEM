#pragma once
#include "local_optimal_scheme.h"
#include "finite_elements_method_intarface.h"
#include "libraries.h"

class one_dimensional_task:finite_elements_method_intarface
{
public:
	one_dimensional_task(std::map<int, int> nodeIndex, std::vector<double> point, std::vector<int> finElemNode[2], std::vector<double> lambda, double boundaryX[2]);
	std::vector<double> runSolver() override;
private:
	void calcLocalMatrix(int cur_el) override; 
	void compileSLAU() override;
	void imposeBoundaryCond() override;

	std::vector<double> _point; // Point coordinate.
	double _boundaryX[2]; // Boundary conditions values
};
