/*************************************************************
	Developed by Alicia Garrido Peña (2020)

	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	
	Please, if you use this implementation cite the two papers above in your work. 
*************************************************************/

#ifndef RAMPGENERATOR_H
#define RAMPGENERATOR_H


/*! RampGenerator class
 * Synapse as described in Vavoulis et al. 
 * Connect two neurons in a gradual synapse.
 */
class RampGenerator
{
		
	double min; ///<Minimun value of the ramp
	double max; ///<Maximum value of the ramp
	double stim_inc; ///<Increment value of the injected current
	double stim_dur; ///<Duration of the same current value in the ramp. 

	public:

		/*! RampGenerator constructor
		* @brief Creates a new ramp generator with all necessary values.
		* @param min Minimun value of the ramp
		* @param max Maximum value of the ramp
		* @param stim_inc Increment value of the injected current
		* @param stim_dur Duration of the same current value in the ramp. 
		*/
		RampGenerator(double min, double max, double stim_inc, double stim_dur);

		/*! RampGenerator empty constructor
		*/
		RampGenerator(){min=0;max=0;stim_inc=0;stim_dur=0;};


		/*!
		* @brief Computes the corresponding value 
		* @param def Default value for the current.
		* @param _time concrete time instant
		* @return Corresponding current value in _time instant. if def != -1, def value is returned.
		*/	
		double get_ext(double def, double _time);


		/*!
		* @brief Prints Ramp Components 
		*/
		void print();

};


#endif
