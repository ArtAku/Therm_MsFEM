#pragma once
#include "libraries.h"

class local_optimal_scheme
{
public:
	local_optimal_scheme(int n1, int elemcount1, type eps1, vector<type> gg1, vector<type> pr1, vector<int> ig1, vector<int> jg1);
	bool iterStep();
	vector<type> getSolve();

	void mulMV(vector<type> vector);
	type mulSKA(vector<type> vector1, vector<type> vector2);
	type calcNorm(vector<type> vector);
private:
	vector<type> _gg;
	vector<int> _ig;
	vector<int> _jg;

	vector<type> _pr;
	vector<type> _r;
	vector<type> _z;
	vector<type> _p;
	type _alpha;
	type _betta;
	type _eps;
	vector<type> _x;
	vector<type> _Ax;

	int _n;
	int _maxiter = 10000;
	int _elemcount;
};

