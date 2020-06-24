/*************************************************************
	Developed by Alicia Garrido Peña (2020)

	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	
	Please, if you use this implementation cite the two papers above in your work. 
*************************************************************/

#include "ramp_generator.h"
#include <iostream>
using namespace std;


RampGenerator::RampGenerator(double min, double max, double stim_inc, double stim_dur)
{
	if(min==-1||max==-1||stim_inc==-1||stim_dur==-1)
		RampGenerator();
	else
	{
		this->min = min;
		this->max = max;
		this->stim_inc = stim_inc;
		this->stim_dur = stim_dur;
		
	}
}


double RampGenerator::get_ext(double def, double _time)
{

	//If there is a default value for the current don't do ramp. 
	if(def !=-1) return def;


	//number of increments necessaries from min to max
	int max_num_inc = (max-min)/stim_inc; 

	//Get the number of stims computed so far.
	int inc_times = _time/stim_dur; 

	//Get the number of stims computed in the current ramp. 
	int inc_times_ramp  = inc_times%max_num_inc;

	double i_ext=0;
	//If ramp is increasing or decreasing:
	//		increasing := number of increments is less than the maximum from min to max in the current ramp.
	//		decreasing := number of increments is greater than the maximum from min to max in the current ramp. 

	if((inc_times/max_num_inc)%2 == 0) //Is increasing
		i_ext = min + inc_times_ramp * stim_inc;

	else if((inc_times/max_num_inc)%2 == 1) //Is decreasing 
		i_ext = max - inc_times_ramp * stim_inc; 

	return i_ext;
	
}


void RampGenerator::print()
{
	cout << "Min "<< min << " Max " << max << " Stim. increment " << stim_inc << " Stim. duration " << stim_dur << endl;
}