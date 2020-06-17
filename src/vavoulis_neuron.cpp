
#include "vavoulis_neuron.h"

#include <iostream>
using namespace std;

VavoulisModel::VavoulisModel(types neu){
	type = neu;

	names[SO]="SO";names[N1M]="N1M";names[N2v]="N2v";names[N3t]="N3t";
	// type, vs, va, p, q,g_ecs,g_eca ,n_syns,conduc_syn,s,r,Esyn,activation_syn, syns[]
	switch(neu)
	{
		case SO: 
			_variables[v] = -65.0; _variables[va] = -65.0; _variables[p] = p_init_SO; _variables[q] = q_init_SO; _variables[h] = h_init; _variables[n] = n_init; 
			params[g_ecs]= g_ec_general; params[g_eca] = g_ec_general;
			
			// names[type]="SO";
			break;


		case N1M:
			_variables[v] = -65.0; _variables[va] = -65.0; _variables[p] = p_init_N1M; _variables[q] = q_init_N1M; _variables[h] = h_init; _variables[n] = n_init; 
			params[g_ecs]= g_ec_general; params[g_eca] = g_ec_general;

			// names[type]="N1M";
			break;

		case N2v:
			_variables[v] = -65.0; _variables[va] = -65.0; _variables[p] = p_init_N2v; _variables[q] = q_init_N2v; _variables[h] = h_init; _variables[n] = n_init; 
			params[g_ecs]= g_ecs_N2v; params[g_eca] = g_eca_N2v;
			
			// names[type]="N2v";

			break;


		case N3t:
			_variables[v] = -65.0; _variables[va] = -65.0; _variables[p] = p_init_N3t; _variables[q] = q_init_N3t; _variables[h] = h_init; _variables[n] = n_init; 
			params[g_ecs]= g_ec_general; params[g_eca] = g_ec_general;
			
			// names[type]="N3t";
			
			break;
		default:
			_variables[v] = -65; _variables[va] = -65; _variables[p] = p_init_SO; _variables[q] = q_init_SO; _variables[h] = h_init; _variables[n] = n_init; 
			params[g_ecs]= g_ec_general; params[g_eca] = g_ec_general;

			// names[type]="None";
			break;


	}	
}
/////////////////////////////////////////////////////////////////
//////////// 			SOMA 		///////////////////////
/////////////////////////////////////////////////////////////////

double VavoulisModel::dp(double _v,double _va,double _p)
{
	double p_inf = 0.0;
	double tau_p = 0.0;
	switch(type)
	{
		case N1M:
			p_inf = 1 / (1+exp((-38.8 - _v)/10));
			tau_p = tau_p_N1M;
			break;

		case N2v:
			p_inf = 1 / (1+exp((-51 - _v)/10.3));
			tau_p = 28.3 + 44.1 * exp(-(((-11.8 - _va)/26.6)*((-11.8 - _va)/26.6))); 
			break;

		case N3t:
			p_inf = 1 / (1+exp((-61.6 - _v)/5.6));
			tau_p = tau_p_N3t;
			break;
		default:
			return 0.0;
	}

	return (p_inf - _p)/tau_p;
}


double VavoulisModel::dq(double _v,double _va,double _q)
{

	double q_inf = 0.0;
	double tau_q =0.0;
	switch(type)
	{
		case N2v:
			q_inf = 1 / (1+exp((-45-_v)/-3));
			tau_q = 187.6 + 637.7 * exp(-(((-9.5-_va)/23.3)*(((-9.5-_va)/23.3))));
			break;
		case N3t:
			q_inf = 1 / (1+exp((-73.2-_v)/-5.1));
			tau_q = tau_q_N3t;
			break;
		default:
			return 0.0;
	}

	return (q_inf - _q)/tau_q;

}


//Channel current based on the current neuron

double VavoulisModel::Ixs(double _v,double _p,double _q)
{

	switch(type)
	{
		case N1M:
			// Iach
			return 200 * _p*_p*_p *(_v+30); 
		case N2v:
			// Inal
			return 2 *_p*_p*_p * _q*(_v-55);
		case N3t:
			// It
			return 3.27 * _p*_p*_p *_q*(_v-80);
		case SO:
			return 0.0;
		default:
			return 0.0;

	}
}



//Differential for Soma Voltage
double VavoulisModel::dvs(double _v,double _va,double _p,double _q, double iext,double isyn)
{

	double ils = _v + 67;
	double ix = Ixs(_v,_p,_q);
	double iecs = params[g_ecs] *(_v-_va);
	
	dv_value = (iext - ils - ix - iecs - isyn)/tau_s;


	return (iext - ils - ix - iecs - isyn)/tau_s;
}


/////////////////////////////////////////////////////////////////
//////////// 			AXON 		///////////////////////
/////////////////////////////////////////////////////////////////

double VavoulisModel::dh(double _va,double _h)
{
	double h_inf = 1/(1 + exp((-55.2 - _va)/-7.1));
	double tau_h = 1.1 + 7.2 * exp(-(pow(((-61.3 - _va)/22.7),2)));
	return (h_inf - _h)/tau_h;

}

double VavoulisModel::dn(double _va,double _n)
{
	double n_inf = 1/(1 + exp((-30 - _va)/17.4));
	double tau_n = 1.1 + 4.6 * exp(-(pow(((-61 - _va)/54.3),2)));
	return (n_inf - _n)/tau_n;

}


double VavoulisModel::Ina (double _va,double _h)
{
	double m = 1/(1+exp((-34.6-_va)/9.6));
	double inaT = 350 * m*m*m * _h * (_va-55);
	return inaT;
}


double VavoulisModel::Ik (double _va,double _n)
{
	double ik = 90 * _n*_n*_n*_n * (_va + 90);
	return ik;
}


//Differential for Axon Voltage
double VavoulisModel::dva(double _v,double _va,double _h,double _n)
{
	double ila = (_va + 67);
	double ina = Ina(_va,_h);
	double ik = Ik(_va,_n);
	double ieca = params[g_eca] * (_va - _v);


	return (-ila - ina - ik - ieca)/tau_a;
}


void VavoulisModel::funcion_simple(double _time, vector<double> &vars, vector<double> &fvec, double iext,double i_syn)
{	
	this->isyn=i_syn;

	fvec[v] = dvs(vars[v],vars[va],vars[p],vars[q],iext,i_syn); 
	fvec[p] = dp(vars[v],vars[va],vars[p]);
	fvec[q] = dq(vars[v],vars[va],vars[q]);
	fvec[va] = dva(vars[v],vars[va],vars[h],vars[n]);
	fvec[h] = dh(vars[va],vars[h]);
	fvec[n] = dn(vars[va],vars[n]);

   return;
}




const char * VavoulisModel::getName()
{
	return names[type];
}

std::vector<double> VavoulisModel::getVariables()
{
	std::vector<double> v(std::begin(_variables), std::end(_variables));
	return v;
}



void VavoulisModel::set_variables(const std::vector<double> &v)
{
	for (int i = 0; i < n_variables && i < (int)v.size(); ++i)
	{
		_variables[i]=v[i];
	}
}




void VavoulisModel::update_variables(double dt, double time, double iext,double i_syn)
{
	std::vector<double>  fvec(n_variables);
	std::vector<double> vars(_variables,_variables+n_variables);

	funcion_simple(time,vars,fvec,iext, i_syn);

	for (int i = 0; i < n_variables; ++i)
	{
		_variables[i]+=fvec[i]*dt;
	}
}

void VavoulisModel::sumValue(int index,double value)
{
	_variables[index]+=value;
}


double VavoulisModel::getVar(int index)
{
	return _variables[index];
}

void VavoulisModel::setVar(int index,double value)
{
	_variables[index] = value;
}

void VavoulisModel::print()
{
	printf("%s %.2f",getName(),V());
}