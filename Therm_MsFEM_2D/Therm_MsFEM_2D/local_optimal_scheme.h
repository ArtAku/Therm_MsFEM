#pragma once
#include "libraries.h"
// Realisation of local optimal scheme - SLAU calculation method. Matrix is stored in sparse line format.


class local_optimal_scheme
{
public:
	local_optimal_scheme(int n1, int elemcount1, double eps1, std::vector<double> gg1, std::vector<double> pr1, std::vector<int> ig1, std::vector<int> jg1);
	std::vector<double> getSolve();
	
private:
	bool iterStep(); // Make one iteration
	void multiplyMV(std::vector<double> vector); // Calculate matrix-vector product. Result in Ax.
	double multiplyScalarVV(std::vector<double> vector1, std::vector<double> vector2); // Calculate scalar product of vectors
	
	int _n;
	int _maxiter = 10000;
	int _elemcount;
	double* _gg;
	int* _ig;
	int* _jg;
	double* _pr;


	double _eps;

	std::vector<double> _r;
	std::vector<double> _z;
	std::vector<double> _p;
	std::vector<double> _x;
	double _alpha;
	double _betta;

	std::vector<double> _MVproduct;
};

