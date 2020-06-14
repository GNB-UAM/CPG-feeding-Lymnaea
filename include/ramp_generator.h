
#ifndef RAMPGENERATOR_H
#define RAMPGENERATOR_H


/*! RampGenerator class
 * Synapse as described in Vavoulis et al. 
 * Connect two neurons in a gradual synapse.
 */
class RampGenerator
{
		
	double min; //<Minimun value of the ramp
	double max; //<Maximum value of the ramp
	double stim_inc; //<Increment value of the injected current
	double stim_dur; //<Duration of the same current value in the ramp. 

	public:

		/*! RampGenerator constructor
		* @brief Creates a new ramp generator with all necessary values.
		* @param min Minimun value of the ramp
		* @param max Maximum value of the ramp
		* @param stim_inc Increment value of the injected current
		* @param stim_dur Duration of the same current value in the ramp. 
		*/
		RampGenerator(double min, double max, double stim_inc, double stim_dur);


		/*!
		* @brief Computes the corresponding value 
		* @param def Default value for the current.
		* @param _time concrete time instant
		* @return Corresponding current value in _time instant. if def != -1, def value is returned.
		*/	
		double get_ext(double def, double _time);

};


#endif
