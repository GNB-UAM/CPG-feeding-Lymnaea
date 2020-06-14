
#ifndef RAMPGENERATOR_H
#define RAMPGENERATOR_H

class RampGenerator
{
		
	double min;
	double max;
	double stim_inc;
	double stim_dur;

	public:
		RampGenerator(double min, double max, double stim_inc, double stim_dur);

		double get_ext(double def, double _time);

};


#endif
