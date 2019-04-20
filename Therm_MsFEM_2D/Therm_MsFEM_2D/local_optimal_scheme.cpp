#include "local_optimal_scheme.h"

local_optimal_scheme::local_optimal_scheme(int n1, int elemcount1, double eps1, std::vector<double> gg1, std::vector<double> pr1, std::vector<int> ig1, std::vector<int> jg1)
{
	_n = n1;
	_elemcount = elemcount1;
	_eps = eps1;
	_x.reserve(_n);
	_p.reserve(_n);
	_r.reserve(_n);
	_z.reserve(_n);
	_MVproduct.reserve(_n);
	_ig = &ig1[0];
	_jg = &jg1[0];
	_gg = &gg1[0];
	_pr = &pr1[0];
	for (int i(0); i < _n; i++)
	{
		_x.push_back(1);
		_r.push_back(0);
		_z.push_back(0);
		_p.push_back(0);
		_MVproduct.push_back(0);
	}
}


bool local_optimal_scheme::iterStep()
{
	int i;
	double S1, S2;
	S1 = multiplyScalarVV(_p, _r);
	S2 = multiplyScalarVV(_p, _p);
	_alpha = S1 / S2;
	for (int i(0); i < _n; i++)
		_x[i] += _alpha * _z[i];
	for (int i(0); i < _n; i++)
		_r[i] -= _alpha * _p[i];
	multiplyMV(_r);
	S1 = multiplyScalarVV(_p, _MVproduct);
	_betta = -(S1 / S2);
	for (int i(0); i < _n; i++)
	{
		_z[i] = _r[i] + _betta * _z[i];
		_p[i] = _MVproduct[i] + _betta * _p[i];
	}
	if (multiplyScalarVV(_r, _r) < _eps)
		return true;
	else return false;
}

std::vector<double> local_optimal_scheme::getSolve()
{
	bool quit = false;
	multiplyMV(_x);
	for (int i(0); i < _n; i++)
	{
		_r[i] = _pr[i] - _MVproduct[i];
		_z[i] = _r[i];
	}
	multiplyMV(_z);
	for (int i(0); i < _n; i++)
		_p[i] = _MVproduct[i];
	for (int i(0); i < _maxiter && !quit; i++)
		quit = iterStep();
	return _x;
}

void local_optimal_scheme::multiplyMV(std::vector<double> vector)
{
	for (int i(0); i < _n; i++)
		_MVproduct[i] = 0;

	for (int i(0); i < _n; i++)
		for (int j = _ig[i], c = _ig[i + 1]; j < c; j++)
			_MVproduct[i] += vector[_jg[j]] * _gg[j];
}


double local_optimal_scheme::multiplyScalarVV(std::vector<double> vector1, std::vector<double> vector2)
{
	int i;
	double res(0);
	for (i = 0; i < _n; i++)
		res += vector1[i] * vector2[i];
	return res;
}