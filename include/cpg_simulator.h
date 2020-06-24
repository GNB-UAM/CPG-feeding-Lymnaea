/*************************************************************
	Developed by Alicia Garrido Peña (2020)

	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	
	Please, if you use this implementation cite the two papers above in your work. 
*************************************************************/


#include "vavoulis_synapse.h"
#include "vavoulis_neuron.h"
#include "ramp_generator.h"

#define SPIKE_TH -50.0
#define MIN_SPIKE_CHANGE 0.001
#define MAX_VARS 18

#define N_NEU 4
#define N_VARS 6

/*! CPGSimulator class
 * Complete circuit class defines neurons and synapses between them to simulate the circuit activity.
 */
class CPGSimulator
{

	std::vector<VavoulisModel > neurons; ///<Vector of neurons 
	int n_neurons; ///<number of neurons initialized
	std::vector<vector< VavoulisSynapse> > syns;///<Vector of synapses vector associated to each neuron
	std::vector<double> c_values; ///<Current value vector (same ids as neurons vector)
	RampGenerator rg; ///<RampGenerator object, contains ramp stimulation
	int connection; ///<Type of connection in the CPG

public:
	/*!Integration methods types
	*/
	enum integrators{EULER,RUNGE,n_integrators};
	CPGSimulator(); ///< Void constructor

	/*! CPGSimulator constructor
	* @brief Creates a new CPG simulator using init method.
	* @param connection type of connection between neurons
	* @param c_values Current value vector (same ids as neurons vector)
	* @param rg RampGenerator object, contains ramp stimulation routines
	*/
	CPGSimulator(int connection, std::vector<double> c_values,RampGenerator rg);

	/*!
	* @brief Assign attributes value depending on the connection type
	* @param connection type of connection between neurons
	* @param c_values Current value vector (same ids as neurons vector)
	* @param rg RampGenerator object, contains ramp stimulation routines
	*/
	int init(int connection, std::vector<double> c_values,RampGenerator rg);

	/*!
	* @brief Simulates activity in the CPG from the initialized CPGSimulator. 
	* @param f File stream
	* @param f_spks Spikes file stream 
	* @param iters Iterations of the simulation
	* @param dt Time step
	* @param integration Integration Method
	* @param satiated_ini Start instant of satiated activity (in iterations)
	* @param satiated_end End instant of satiated activity (in iterations)
	*/
	void simulate(FILE * f,FILE * f_spks,double iters,double dt,integrators integration,double statiated_ini,double satiated_end);

	/*!
	* @brief Prints CPG components: All neurons and synapses initialized
	*/
	void print();


private:
	
	/*!
	* @brief General update function, this function call either update_euler or update_runge
	* @param _time Current time instant
	* @param integr Integration Method
	* @param dt Time step
	*/
	double update_all(double _time, integrators integr,double dt);
	/*!
	* 	@brief Euler update function, updates
	* @param _time Current time instant
	* @param dt Time step
	*/
	void update_euler(double _time, double dt);
	/*!
	* @brief Runge-Kutta update function, updates function using intey and diffs auxiliary function
	* @param _time Current time instant
	* @param dt Time step
	*/
	void update_runge(double _time, double dt);
	/*!
	* 	@brief Intey auxiliar function, obtains general vectors with the result of each differential equation for each neuron variables and its associated synapses.
	*	@param _time Current time instant
	* 	@param v_variables general vector matrix with all variables (neuron+synapses)
	* 	@param v_fvec general return vector matrix with all differential equations value for each variable (neuron+synapses)
	*/
	void diffs(double _time,const std::vector<std::vector<double > > & v_variables, std::vector<std::vector<double > >  &v_fvec);
	/*!
	* 	@brief Performs Runge-Kutta integration with middle steps. New variable defined with a "global" vector containing both neurons and synapses.
	* 	Complete array example with N1M-N2v-N3t connection:
	* 		SO neuron variables (no synapses):					 v,va,p,q,h,n
	* 		N1M neuron variables and synapses N1M-N2v;N1M-N3t:	 v,va,p,q,h,n,s,r,s,r
	* 		N2v variables and synapses N2v-N1M:					 v,va,p,q,h,n,s,r
	* 		N3t variables and synapses N3t-N1M;N3t-N2v:			 v,va,p,q,h,n,s,r,s,r
	*
	*	@param _time Current time instant
	*	@param inc_integracion Time step
	*/
	double intey(double _time, double inc_integracion);
	/*!
	* @brief Detect possible spikes in each neuron and writes it in the associated spike file. When no spike is found ',' is written in the corresponding column.
	* @param f_spks Spikes file stream 
	* @param prevs vector of derivate previous values 
	* @param t time value 
	*/
	void detect_spikes(FILE * f_spks, std::vector<double> &prevs, double t );

	/*!
	* @brief Writes V value and Isyn or injected current depending on syns_flag value
	* @param f File stream
	* @param t time instant
	* @param c current value
	*/
	void write(FILE *f,double t,double c);

};