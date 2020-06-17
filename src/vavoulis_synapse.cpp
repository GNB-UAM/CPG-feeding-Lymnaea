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



void VavoulisSynapse::funcion_simple(double time, const std::vector<double> & vars, std::vector<double> &fvec, double vpre)
{
	fvec[s] = ds(vars[r],vars[s]); 
	fvec[r] = dr(vpre,vars[r]);

   return;
}

void VavoulisSynapse::funcion_simple(double time, std::vector<double> &fvec, double vpre)
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

	funcion_simple(time,vars,fvec,vpre);

	for (int i = 0; i < n_variables; ++i)
	{
		_variables[i]+=fvec[i]*dt;
	}
}




void VavoulisSynapse::update_variables(double dt, double time)
{
	update_variables(dt,time,pre->V());
}


void VavoulisSynapse::showParams()
{
	for (int i=0; i<n_params; i++)
		cout <<params[i] << " ";
}

void VavoulisSynapse::print()
{
	cout << "Pos: ";
	pos->print(); cout<< endl;
	cout << "Pre: "; pre->print();  cout<< endl;
}



// void VavoulisSynapse::update(VavoulisModel::integrators integr,double dt)
// {
// 	if(integr == VavoulisModel::EULER)
// 	{
// 		update_euler(dt);
// 	}
// 	else if(integr == VavoulisModel::RUNGE)
// 	{

// 		update_runge(dt);
// 	}
// }


// void VavoulisSynapse::update_runge(double dt)
// {
// 	intey(dt);

// }


// void VavoulisSynapse::update_euler(double dt)
// {
// 	double aux_variables[n_variables];

// 	aux_variables[s] = _variables[s] + dt*ds(_variables);
// 	aux_variables[r] = _variables[r] + dt*dr(_variables);

// 	std::copy(aux_variables,aux_variables+n_variables,_variables);

// }



// void VavoulisSynapse::funcion(double * variables, double * fvec)
// {

// 	fvec[s] = ds(); 
// 	fvec[r] = dr();

//    return;
// }



// void VavoulisSynapse::print_params()
// {
// 	for(int i=0; i < n_params; i++)
// 	{
// 		cout << params[i] << endl;
// 	}
// }


// double VavoulisSynapse::dr(double * vars)
// {
// 	//Steady value: Vpre sigmoid function
// 	// double r_inf = 1/(1 + exp((RUMBRAL - Vpre)/R_TAU));
// 	double r_inf = 1/(1 + exp((-40 - pre->V())/2.5));
// 	double r_value = (r_inf - vars[r])/params[activation_syn];
// 	return r_value;
// }


// double VavoulisSynapse::ds(double * vars)
// {
// 	return (vars[r] - vars[s]) / params[activation_syn];
// }


// double VavoulisSynapse::Isyn()
// {
// 	return params[conduc_syn] * _variables[s] * (pos->V() - params[Esyn]);
// }

