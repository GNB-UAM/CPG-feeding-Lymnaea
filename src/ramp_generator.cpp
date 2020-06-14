
#include "ramp_generator.h"
#include <iostream>
using namespace std;


	RampGenerator::RampGenerator(double min, double max, double stim_inc, double stim_dur)
	{
		this->min = min;
		this->max = max;
		this->stim_inc = stim_inc;
		this->stim_dur = stim_dur;
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

		// cout << _time << " "<< i_ext << endl;

		return i_ext;
		
	}
