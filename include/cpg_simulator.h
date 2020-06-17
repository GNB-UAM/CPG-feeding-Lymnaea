#include "vavoulis_synapse.h"
#include "vavoulis_neuron.h"
#include "ramp_generator.h"

#define SPIKE_TH -50.0
#define MIN_SPIKE_CHANGE 0.001
#define MAX_VARS 18

#define N_NEU 4
#define N_VARS 6
class CPGSimulator
{

	std::vector<VavoulisModel > neurons;
	int n_neurons;
	std::vector<vector< VavoulisSynapse> > syns;
	std::vector<double> c_values;


public:
	enum integrators{EULER,RUNGE,n_integrators};
	CPGSimulator();
	CPGSimulator(int connection, std::vector<double> c_values);

	int init(int connection, std::vector<double> c_values);
	void simulate(FILE * f,FILE * f_spks,double iters,double dt,integrators integration,double feed_ini,double feed_end,int connection);
	void print();


private:
	RampGenerator rg;
	
	double update_all(double _time, integrators integr,double dt);
	double update_euler(double _time, double dt);
	void update_runge(double _time, double dt);
	void funcion(double _time,const std::vector<std::vector<double > > & v_variables, std::vector<std::vector<double > >  &v_fvec);
	double intey(double _time, double inc_integracion);
	void detect_spikes(FILE * f_spks, std::vector<double> &prevs, double t );

	void write(FILE *f,double t,double c,int syns_flag);

};