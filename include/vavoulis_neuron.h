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
	enum vars_names {v,va,p,q,h,n,n_variables}; //< Varibles names references for _variables array
	enum pars_names {g_ecs,g_eca,n_params}; //< Parameters names references for _params array
	double _variables[n_variables]; //< Variables array. All elements in this array have a differential equation associated in the form of dvar(...)
	double params[n_params]; //< Parameter array
	double dv_value; //< Last dv_value computed
	double isyn; //< Last synaptic input received
	double n_syns; //< Number of synapses received by this neuron. 

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
	// enum integrators{EULER,RUNGE,n_integrators};
	types type; ///< Neuron type from types 

	// std::vector<void *> synapses;

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

	// void update(integrators integr,double dt,double i_ext,double isyn);


	// void add_synapse(VavoulisSynapse * synapse);
	// void update_synapses(integrators integr,double dt);
	// double Isyn();

	// int n_syns;
	const char * getName();

	// void addSynapse(void * syn);
	void sumValue(int index,double value);
	double getVar(int index); 
	void setVar(int index,double value);


	double getNSynapses(){return n_syns;}//< Number of synapses getter. 
	// double getNVarsTotal();
	// double getNVars(){return n_variables;}
	
	static int getNVars(){return n_variables;}//< returns the number of variables. 

	/*!
	 * 
	 * @brief variables array getter
	 * @return A copy of variables array as a std vector.
	 * @see Vavoulis et al. 
	 */
	std::vector<double> getVariables();

	/*!
	 * 
	 * @brief Updates the n_variables in variables array with the data in v. 
	 * @param v vector containing new variables values.
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. 
	 */
	void set_variables(const std::vector<double> &v);

	// void funcion_simple(double _time, double *vars, double *fvec,double v_vars[n_types] , double iext,RampGenerator * rg);
	// void funcion_simple(double _time, double *vars, double *fvec,double v_vars[n_types][18] , double iext,RampGenerator * rg);
	// void update_variables(double dt, double time, double iext,double i_syn);
	void update_variables(double dt, double time, double iext,double i_syn);
	/*!
	 * 
	 * @brief Returns an array with all values obtained by the differential equations
	 * @param _time Current time for the differential equation. 
	 * @param vars vector with previous instant variables vales. 
	 * @param fvec return vector with computed differential values.
	 * @param iext injected current received.
	 * @param isyn synaptic current received.
	 * @param rg RampGenerator object. 
	 * @see Vavoulis et al. 
	 */
	void funcion_simple(double _time, vector<double> &vars, vector<double> &fvec, double iext,double i_syn);
	
	// double update_euler(double _time, double dt, double iext,RampGenerator * rg);

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
	 * @see Vavoulis et al. 
	 */
	double dp(double _v,double _va,double _p); 

	/*!
	 * q differential equation
	 * @brief differential equation corresponding to q neurotransmitor in iX channels from the somatic compartment. Its equation depends on the neuron type. 
	 * @param _v Voltage value at soma compartment
	 * @param _va Voltage value at axon compartment
	 * @param _q previous q value. (usually _variables[q])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. 
	 */
	double dq(double _v,double _va,double _q);


	/*!
	 * Slow channel from the somatic compartment
	 * @brief Ixs channel represents either IACh, INaL or IT for N1M, N2v or N3t, respectively. Represent slow activity in the somatic compartment. 
	 * @param _v Voltage value at soma compartment
	 * @param _p value (usually _variables[p])
	 * @param _q value (usually _variables[q])
	 * @return channel value in mV. 
	 * @see Vavoulis et al. 
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
	 * @see Vavoulis et al. 
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
	 * @see Vavoulis et al. 
	 */
	double dh(double _va,double _h);

	/*!
	 * n differential equation
	 * @brief differential equation corresponding to n neurotransmitor in iK from the axonal compartment. 
	 * @param _va Voltage value at axon compartment
	 * @param _n previous n value. (usually _variables[n])
	 * @return result for the equation for the corresponding neuron. 
	 * @see Vavoulis et al. 
	 */
	double dn(double _va,double _n);
	
	/*!
	 * Fast channel INaT from the axonal compartment
	 * @brief INaT channel represents an inactivating sodium current, part of fast axonal compartment.
	 * @param _va Voltage value at axonal compartment
	 * @param _h value (usually _variables[h])
	 * @return channel value in mV. 
	 * @see Vavoulis et al. 
	 */
	double Ina (double _va,double _h);
	/*!
	 * Fast channel IK from the axonal compartment
	 * @brief IK channel represents an rectifier potassium current, part of fast axonal compartment.
	 * @param _va Voltage value at axonal compartment.
	 * @param _n value (usually _variables[n])
	 * @return channel value in mV. 
	 * @see Vavoulis et al. 
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
	 * @see Vavoulis et al. 
	 */
	double dva(double _v,double _va,double _h,double _n);


	// void funcion(double * variables, double * fvec, int iext,double isyn);
	// double intey(double inc_integracion, double iext, double isyn);

	// double update_runge(double dt, double iext,double isyn);

	// double Isyn_total(double *);



};

#endif 