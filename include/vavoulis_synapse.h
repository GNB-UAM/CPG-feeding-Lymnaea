#ifndef SYNS_H
#define SYNS_H


#define r_init 0.000045398
#define s_init 0.000045398

#define EXCIT 0.0
#define INHIB -90.0

#define SLOW 200.0
#define FAST 50.0

#include"vavoulis_neuron.h"



/*! VavoulisSynapse class
 * Synapse as described in Vavoulis et al. 
 * Connect two neurons in a gradual synapse.
 */
class VavoulisSynapse{

	enum vars_names {s,r,n_variables}; ///< Varibles names references for _variables array
	enum params{conduc_syn, activation_syn,Esyn,n_params}; ///< Parameters names references for _params array

	double _variables[n_variables]; ///< Variables array. All elements in this array have a differential equation associated in the form of dvar(...)
	double params[n_params];  ///< Parameter array
	int pre_type; ///< Presynaptic neuron type.

	VavoulisModel *pos; ///< posynaptic neuron reference
	VavoulisModel *pre;///< pre-synaptic neuron reference

public:

	/*! VavoulisSynapse constructor
	* @brief Creates a new synapse model with all the corresponding params and variables initial values assigned.
	* @param n1 Pointer to pos-synaptic neuron
	* @param n2 Pointer to pre-synaptic neuron
	* @param _conduc_syn Maximal synaptic conductance parameter value
	* @param _act_syn Activation time constant parameter value
	* @param _Esyn reversal potential parameter value
	*/
	VavoulisSynapse(VavoulisModel *n1, VavoulisModel* n2,double _conduc_syn,double _act_syn,double _Esyn);

	/*! VavoulisSynapse void constructor
	* @brief Creates a new synapse model with all the corresponding params and variables initial values initialized to 0.
	*/
	VavoulisSynapse(); 
	

	static int getNVars(){return n_variables;} ///< returns the number of variables. 
	int getPreType() {return pre_type;} ///< Presynaptic neuron type getter
	
	/*!
	 * 
	 * @brief variables array getter
	 * @return A copy of variables array as a std vector.
	 * @see Vavoulis et al. 
	 */
	std::vector<double> getVariables();


	/*! 
	* @brief Sets params values of the synapse.
	* @param n1 Pointer to pos-synaptic neuron
	* @param n2 Pointer to pre-synaptic neuron
	* @param _conduc_syn Maximal synaptic conductance parameter value
	* @param _act_syn Activation time constant parameter value
	* @param _Esyn reversal potential parameter value
	*/
	void setSynapse(VavoulisModel *n1, VavoulisModel* n2,double _conduc_syn,double _act_syn,double _Esyn);

	/*!
	 * 
	 * @brief Updates the n_variables in variables array with the data in v. 
	 * @param v vector containing new variables values.
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. 
	 */
	void set_variables(const std::vector<double> &v);

	/*!
	 * 
	 * @brief Returns an array with all values obtained by the differential equations
	 * @param time Current time for the differential equation. 
	 * @param vars vector with previous instant variables vales. 
	 * @param fvec return vector with computed differential values.
	 * @param vpre Voltage vale in the somatic compartment from the Presynaptic neuron.
	 * @see Vavoulis et al. 
	 */
	void funcion_simple( double time, const std::vector<double> & vars, std::vector<double> &fvec, double vpre);

	void funcion_simple(double time, std::vector<double> &fvec, double vpre);

	void update_variables(double dt, double time, double vpre);
	void update_variables(double dt, double time);



	/*!
	 * Synaptic current value
	 * @brief Synaptic current resulted from Pos- and presynaptic neuron.
	 * @param _v Soma voltage value from pos-synaptic neuron.
	 * @param _r previous r value. (usually _variables[r])
	 * @return resulting synaptic value in mV.
	 * @see Vavoulis et al. 
	 */
	double Isyn(const vector<double> &variables,double v);
	double Isyn();
	double Isyn(double v);
	// double Isyn_runge(double _s,double _v);
	void print_params();



	void showParams();
	void update_euler(double dt);

private:

	// double dr(double * );
	// double ds(double * );

	/*!
	 * r differential equation
	 * @brief differential equation corresponding to r steady state in the synaptic equation.
	 * @param _vpre Voltage vale in the somatic compartment from the Presynaptic neuron.
	 * @param _r previous r value. (usually _variables[r])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. 
	 */
	double dr(double _vpre,double _r);

	/*!
	 * s differential equation
	 * @brief differential equation corresponding to s activation value in the synaptic equation.
	 * @param _vpre Voltage vale in the somatic compartment from the Presynaptic neuron.
	 * @param _r value of r. (usually _variables[r])
	 * @param _s previous s value. (usually _variables[s])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. 
	 */
	double ds(double _r,double _s);
	// double Isyn();

};
#endif 