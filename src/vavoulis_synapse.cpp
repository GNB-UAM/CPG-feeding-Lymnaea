/*************************************************************
	Developed by Alicia Garrido Peña (2020)

	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	
	Please, if you use this implementation cite the two papers above in your work. 
*************************************************************/


#include "vavoulis_synapse.h"
#include <iostream>
using namespace std;

VavoulisSynapse::VavoulisSynapse(VavoulisModel *n1, VavoulisModel* n2,double _conduc_syn,double _activation_syn,double _Esyn)
{
	pos = n1;
	pre = n2;
	params[conduc_syn] = _conduc_syn;
	params[activation_syn] = _activation_syn;
	params[Esyn] = _Esyn;

	_variables[s] = s_init;
	_variables[r] = r_init;

	pre_type = pre->type;
}


VavoulisSynapse::VavoulisSynapse()
{
	pos = 0;
	pre = 0;
	params[conduc_syn] = 0;
	params[activation_syn] = 0;
	params[Esyn] = 0;

	_variables[s] = 0;
	_variables[r] = 0;
	pre_type=0;
}


void VavoulisSynapse::setSynapse(VavoulisModel *n1, VavoulisModel* n2,double _conduc_syn,double _act_syn,double _Esyn)
{
	pos = n1;
	pre = n2;
	params[conduc_syn] = _conduc_syn;
	params[activation_syn] = _act_syn;
	params[Esyn] = _Esyn;

	_variables[s] = s_init;
	_variables[r] = r_init;
	pre_type = pre->type;
}

double VavoulisSynapse::dr(double _vpre,double _r)
{
	double r_inf = 1/(1 + exp((-40 - _vpre)/2.5));
	double r_value = (r_inf - _r)/params[activation_syn];
	return r_value;
}


double VavoulisSynapse::ds(double _r,double _s)
{
	return (_r-_s) / params[activation_syn];
}

double VavoulisSynapse::Isyn(const vector<double> &variables,double v)
{
	return params[conduc_syn] * variables[s] * (v - params[Esyn]);
}

double VavoulisSynapse::Isyn(double v)
{
	return params[conduc_syn] * _variables[s] * (v - params[Esyn]);
}


double VavoulisSynapse::Isyn()
{
	return params[conduc_syn] * _variables[s] * (pos->V() - params[Esyn]);
}



void VavoulisSynapse::diffs_fun(double time, const std::vector<double> & vars, std::vector<double> &fvec, double vpre)
{
	fvec[s] = ds(vars[r],vars[s]); 
	fvec[r] = dr(vpre,vars[r]);

   return;
}

void VavoulisSynapse::diffs_fun(double time, std::vector<double> &fvec, double vpre)
{
	fvec[s] = ds(_variables[r],_variables[s]); 
	fvec[r] = dr(vpre,_variables[r]);

   return;
}

std::vector<double> VavoulisSynapse::getVariables()
{
	std::vector<double> v(std::begin(_variables), std::end(_variables));
	return v;
}

void VavoulisSynapse::set_variables(const std::vector<double>& v)
{
	// cout <<"N_VARS"<< " "<<  n_variables<<endl;
	for (int i = 0; i < n_variables; ++i)
	{
		_variables[i]=v[i];
	}
}

void VavoulisSynapse::update_variables(double dt, double time, double vpre)
{
	std::vector<double>  fvec(n_variables);
	std::vector<double> vars(_variables,_variables+n_variables);

	diffs_fun(time,vars,fvec,vpre);

	for (int i = 0; i < n_variables; ++i)
	{
		_variables[i]+=fvec[i]*dt;
	}
}




void VavoulisSynapse::update_variables(double dt, double time)
{
	update_variables(dt,time,pre->V());
}


void VavoulisSynapse::print_params()
{
	for (int i=0; i<n_params; i++)
		cout <<params[i] << " ";
}

void VavoulisSynapse::print()
{
	cout << "\tPos: ";
	pos->print(); cout<< endl;
	cout << "\tPre: "; pre->print();  cout<< endl;
}


