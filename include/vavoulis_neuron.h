/*************************************************************
	Developed by Alicia Garrido Peña (2020)

	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	
	Please, if you use this implementation cite the two papers above in your work. 
*************************************************************/


#include <vector>
#include <math.h>
#include <iterator>

using namespace std;

#ifndef MODEL_H
#define MODEL_H

//All values bellow are set out in Tables 1,2 and 3 in Vavoulis et al.
//Init values are obtained as the result of their corresponding differential equation.

//------------------Model parameters-------------------------------
#define tau_p_N1M 250
#define tau_q_N1M 0
#define tau_p_N2v 0
#define tau_q_N2v 0
#define tau_p_N3t 4
#define tau_q_N3t 400
#define tau_p_SO 0
#define tau_q_SO 0

#define tau_s 10
#define tau_a 10


#define g_eca_N2v 0.06
#define g_ecs_N2v 0.55
#define g_ec_general 8

#define h_init 0.799
#define n_init 0.118

#define h_init_30 0.000006144
#define n_init_30 0.96917

#define p_init_N1M 0.0678
#define q_init_N1M 0 //not necessary

#define p_init_SO 0 //not necessary
#define q_init_SO 0 //not necessary

#define p_init_N2v 0.2043
#define q_init_N2v 0.3527

#define p_init_N3t 0.3527
#define q_init_N3t 0.1668
//-------------------------------------------------------------------


/*! VavoulisModel class
 * Single neuron class following Vavoulis et al. description.
 */
class VavoulisModel{
	/*!
	 Varibles names references for _variables array
	*/
	enum vars_names {v,va,p,q,h,n,n_variables}; 
	/*!
	 Parameters names references for _params array
	*/
	enum pars_names {g_ecs,g_eca,n_params}; 
	double _variables[n_variables]; ///< Variables array. All elements in this array have a differential equation associated in the form of dvar(...)
	double params[n_params]; ///< Parameter array
	double dv_value; ///< Last dv_value computed
	double isyn; ///< Last synaptic input received
	double n_syns; ///< Number of synapses received by this neuron. 

public:
	/*!
 	* @brief Enumeration type with the 4 possible neuron types of this model. 
 	* 	SO -> Slow oscillator modulatort interneuron. 
 	* 	N3t -> N3 tonic.
 	* 	N2v -> N2 ventral.
 	* 	N1M -> N1 medial.
 	*/
	enum types{SO,N1M,N2v,N3t,n_types};

	const char * names[n_types]; ///< Array with each neuron name. 
	types type; ///< Neuron type from types 


	/*! VavoulisModel constructor
	* @brief Creates a new neuron model with all the corresponding params and variables initial values depending on its type.
	* @param neu neuron type from enum types.
	*/
	VavoulisModel(types neu);


	/*!
 	* @brief Voltage soma value getter
 	* @return v value in _variables
 	*/
	double V(){return _variables[v];}
	/*!
 	* @brief Voltage axon value getter
 	* @return va value in _variables
 	*/
	double Va(){return _variables[va];}

	/*!
 	* @brief V derivative getter.
 	* @return Last result of dv(). 
 	* @see dvs(double _v,double _va,double _p,double _q, double iext,double isyn);
 	*/
	double getdV(){return dv_value;}

	/*!
 	* @brief Isyn getter.
 	* @return Last value of isyn received. 
 	* @see dvs
 	*/
	double getIsyn(){return isyn;}

	
	/*!
 	* @brief Name getter
 	* @return Neuron name string 
 	*/
	const char * getName();



	double getNSynapses(){return n_syns;}///< Number of synapses getter. 
	
	static int getNVars(){return n_variables;}///< returns the number of variables. 

	void sumValue(int index,double value);///< Sets _variables[index]+=value
	double getVar(int index); ///< Returns _variables[index]
	void setVar(int index,double value);///< Sets _variables[index]=value
	
	/*!
	 * 
	 * @brief variables array getter
	 * @return A copy of variables array as a std vector.
	 */
	std::vector<double> getVariables();

	/*!
	 * 
	 * @brief Updates the n_variables in variables array with the data in v. 
	 * @param v vector containing new variables values.
	 * @return result for the equation for the corresponding neuron. 
	 */
	void set_variables(const std::vector<double> &v);

	
	/*!
	 * @brief Integrates variables based on dt.
	 * @param dt Integration time step
	 * @param time Current time for the differential equation. 
	 * @param iext injected current received.
	 * @param i_syn synaptic current received.
	 */
	void update_variables(double dt, double time, double iext,double i_syn);

	/*!
	 * 
	 * @brief Returns an array with all values obtained by the differential equations
	 * @param _time Current time for the differential equation. 
	 * @param vars vector with previous instant variables values. 
	 * @param fvec return vector with computed differential values.
	 * @param iext injected current received.
	 * @param i_syn synaptic current received.
	 */
	void diffs_fun(double _time, vector<double> &vars, vector<double> &fvec, double iext,double i_syn);
	

	/*!
	*
	* @brief Prints Neuron.
	* Format: 
	*  		Neuron-Name V-value
	*/
	void print();
private: 
		/////////////////////////////////////////////////////////////////
	//////////// 			SOMA 		///////////////////////
	/////////////////////////////////////////////////////////////////

	/*!
	 * p differential equation
	 * @brief differential equation corresponding to p neurotransmitor in iX channels from the Somatic compartment. Its equation depends on the neuron type. 
	 * @param _v Voltage value at soma compartment
	 * @param _va Voltage value at axon compartment
	 * @param _p previous p value. (usually _variables[p])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. [1] 
	 */
	double dp(double _v,double _va,double _p); 

	/*!
	 * q differential equation
	 * @brief differential equation corresponding to q neurotransmitor in iX channels from the somatic compartment. Its equation depends on the neuron type. 
	 * @param _v Voltage value at soma compartment
	 * @param _va Voltage value at axon compartment
	 * @param _q previous q value. (usually _variables[q])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. [1] 
	 */
	double dq(double _v,double _va,double _q);


	/*!
	 * Slow channel from the somatic compartment
	 * @brief Ixs channel represents either IACh, INaL or IT for N1M, N2v or N3t, respectively. Represent slow activity in the somatic compartment. 
	 * @param _v Voltage value at soma compartment
	 * @param _p value (usually _variables[p])
	 * @param _q value (usually _variables[q])
	 * @return channel value in mV. 
	 * @see Vavoulis et al. [1] 
	 */
	//Channel current based on the current neuron
	double Ixs(double _v,double _p,double _q);

	/*!
	 * Differential equation for somatic compartment voltage where slow activity properties are hosted 
	 * @brief Differential equation for voltage in soma. 
	 * @param _v Voltage value at soma compartment.
	 * @param _va Voltage value at axon compartment.
	 * @param _p value (usually _variables[p]).
	 * @param _q value (usually _variables[q]).
	 * @param iext injected current received.
	 * @param isyn synaptic current received.
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. [1] 
	 */
	double dvs(double _v,double _va,double _p,double _q, double iext,double isyn);

	/////////////////////////////////////////////////////////////////
	//////////// 			AXON 		///////////////////////
	/////////////////////////////////////////////////////////////////

	/*!
	 * h differential equation
	 * @brief differential equation corresponding to h neurotransmitor in iNaT from the axial compartment. 
	 * @param _va Voltage value at axon compartment
	 * @param _h previous h value. (usually _variables[h])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. [1] 
	 */
	double dh(double _va,double _h);

	/*!
	 * n differential equation
	 * @brief differential equation corresponding to n neurotransmitor in iK from the axonal compartment. 
	 * @param _va Voltage value at axon compartment
	 * @param _n previous n value. (usually _variables[n])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. [1] 
	 */
	double dn(double _va,double _n);
	
	/*!
	 * Fast channel INaT from the axonal compartment
	 * @brief INaT channel represents an inactivating sodium current, part of fast axonal compartment.
	 * @param _va Voltage value at axonal compartment
	 * @param _h value (usually _variables[h])
	 * @return channel value in mV. 
	 * @see Vavoulis et al. [1] 
	 */
	double Ina (double _va,double _h);
	/*!
	 * Fast channel IK from the axonal compartment
	 * @brief IK channel represents an rectifier potassium current, part of fast axonal compartment.
	 * @param _va Voltage value at axonal compartment.
	 * @param _n value (usually _variables[n])
	 * @return channel value in mV. 
	 * @see Vavoulis et al. [1] 
	 */
	double Ik (double _va,double _n);

	/*!Differential equation for voltage in axon. 
	 *  
	 * @brief Differential equation for axonal compartment voltage where fast activity properties are hosted.
	 *  Spiking frequency depends on this compartment. Conductance in electric coupling with soma influences the spiking. 
	 * @param _v Voltage value at soma compartment.
	 * @param _va Voltage value at axon compartment.
	 * @param _h value (usually _variables[h]).
	 * @param _n value (usually _variables[n]).
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. [1] 
	 */
	double dva(double _v,double _va,double _h,double _n);


};

#endif 