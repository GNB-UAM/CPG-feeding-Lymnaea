/*
	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	Please, if you use this implementation cite the two papers above in your work. 
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <string>
#include <vector>
#include <iostream>

#include "ramp_generator.h"
#include "vavoulis_synapse.h"

using namespace std;

#define n_variables 8

#define MAX_VARS 18

#define N_NEU 4
#define N_VARS 6


#define SPIKE_TH -50.0
#define MIN_SPIKE_CHANGE 0.001



#define MAX_STRING 20000

#define ERROR 0
#define OK 1

enum prim_types{String,Integer, Float, Double, IntegrationMeth};//<Data types for parsing function

string arg_names[]={"-connection","-file_name","-integrator","-dt","-c_so","-c_n1m","-c_n2v","-c_n3t","-stim_dur","-stim_inc","-MIN_c","-MAX_c","-secs_dur","-rounds","-feed_ini","-feed_end"};//<Arguments possible names
prim_types arg_types[]={Integer,String,IntegrationMeth,Double,Double,Double,Double,Double,Double,Double,Double,Double,Double,Integer,Double,Double}; //Arguments corresponding types

string format = "Format: ./lymn -connection val -file_name val -integration_method -dt val val -c_so val -c_n1m val -c_n2v val -c_n3t val -stim_dur val -stim_inc val -MIN_c val -MAX_c val [-secs_dur val] -rounds val -feed_ini val -feed_end val\n";//<Input format


string methods[] = {"Euler","Runge-Kutta"}; //<Integrator names in String
string headers[] = {"t SO N1M N2v N3t c", "t N1M N2v","t N1M N2v N3t","t SO N1M N2v N3t c","t SO IsynSO N1M IsynN1M N2v IsynN2v N3t IsynN3t",
		"t SO N1M N2v N3t", "t SO N1M N2v N3t c"};//<File headers depending on the connection. 


enum integrators{EULER,RUNGE,n_integrators}; //<Integration types available

/*!
* @brief Parse input arguments defined by its types above in arg_names and arg_types. 
*
*/
int parse_input(int argc,char *argv[], void ** arguments);
/*!
* Prompt help with parameters description
*/
void show_help();

/*!
* 	Detect possible spikes in each neuron and writes it in the associated spike file. When no spike is found ',' is written in the corresponding column.
* @param f_spks Spikes file stream 
* @param neurons vector of neurons 
* @param prevs vector of derivate previous values 
* @param t time value 
*/
void detect_spikes(FILE * f_spks, std::vector<VavoulisModel *> const & neurons, std::vector<double> &prevs, double t );

/*!
* 	General update function, this function call either update_euler or update_runge
* @param _time Current time instant
* @param neurons Vector of neurons 
* @param syns Vector of synapses vector associated to each neuron
* @param integr Integration Method
* @param dt Time step
* @param i_ext Current value vector (same ids as neurons vector)
* @param rg RampGenerator object, contains ramp stimulation
* @return Injected current value of the neuron stimulated (-1)
*/
double update_all(double _time, const vector<VavoulisModel *>& neurons,const std::vector<vector<VavoulisSynapse * > >& syns, integrators integr,double dt,const std::vector<double> & i_ext, RampGenerator * rg);

/*!
* 	@brief Euler update function, updates
* 	@param see update all parameters description
*/
void update_euler(double _time, double dt,std::vector<VavoulisModel *> neurons, std::vector<std::vector<VavoulisSynapse* > > syns ,const std::vector<double> &iext,RampGenerator * rg);
/*!
* 	@brief Runge-Kutta update function, updates function using intey and funcion auxiliar function
* 	@param see update all parameters description
*/
void update_runge(double _time,const std::vector<VavoulisModel * >& neurons,const std::vector<vector<VavoulisSynapse * > >& syns, double dt , const std::vector<double> &iext, RampGenerator * rg);
/*!
* 	@brief Performs Runge-Kutta integration with middle steps. New variable defined with a "global" vector containing both neurons and synapses.
* 	Complete array example with N1M-N2v-N3t connection:
* 		SO neuron variables (no synapses):					 v,va,p,q,h,n
* 		N1M neuron variables and synapses N1M-N2v;N1M-N3t:	 v,va,p,q,h,n,s,r,s,r
* 		N2v variables and synapses N2v-N1M:					 v,va,p,q,h,n,s,r
* 		N3t variables and synapses N3t-N1M;N3t-N2v:			 v,va,p,q,h,n,s,r,s,r
*
* 	@param see update all parameters description
*/
double intey(double _time, const std::vector<VavoulisModel *> &neurons,const std::vector<vector<VavoulisSynapse *> >& syns, double inc_integracion, const std::vector<double> & iext, RampGenerator * rg);
/*!
* 	@brief Intey auxiliar function, obtains general vectors with the result of each differential equation for each neuron variables and its associated synapses.
*	@param _time Current time instant
* 	@param v_variables general vector matrix with all variables (neuron+synapses)
* 	@param v_fvec general return vector matrix with all differential equations value for each variable (neuron+synapses)
*	@param neurons Vector of neurons 
*	@param syns Vector of synapses vector associated to each neuron
*	@param dt Time step
*	@param i_ext Current value vector (same ids as neurons vector)
* 	@param rg RampGenerator object, contains ramp stimulation
* 	@param see update all parameters description
*/
void funcion(double _time,const std::vector<std::vector<double > > & v_variables, std::vector<std::vector<double > >  &v_fvec, const std::vector<VavoulisModel *>& neurons,const std::vector<vector<VavoulisSynapse *> >& syns , const std::vector<double> & iext, RampGenerator * rg);



int main(int argc, char * argv[])
{

	int connection =3;
	char *file_name;
	char file_spikes[MAX_STRING];
	char file_ext[MAX_STRING];
	double secs_dur = -1;
	int iters =-1;
	integrators integration = EULER;
	double dt=0.01;
	double c_so=-1,c_n1m=-1,c_n2v=-1,c_n3t=-1;
	double stim_dur=-1;
	double stim_inc=-1,MIN_c=-1,MAX_c=-1;
	double feed_ini =  -1,feed_end =  -1;

	int rounds=4;

	FILE *f,*f_spks;

	const char * header;

	///////////////////////////////////////
	//Input Parameters 
	///////////////////////////////////////

	if(argc == 1){
		cout <<format<< endl;
		cout << "If any of those parameters is skiped, the default value will be assigned\n ./lymn --help for more info" << endl;
		return -1;
	}
	else if(argc==2 && strcmp(argv[1],"--help")==0)
	{
		show_help();
		return -1;
	}
	else
	{
		void * arguments[] = {&connection,&file_name,&integration,&dt,&c_so,&c_n1m,&c_n2v,&c_n3t,&stim_dur,&stim_inc,&MIN_c,&MAX_c,&secs_dur,&rounds,&feed_ini,&feed_end};
		if(parse_input(argc,argv,arguments)==ERROR)
		{
			cerr << "Error parsing input"<< endl;
			return -1;
		}
		if(!file_name)
		{
			cerr<< "No file name specified"<<endl;
			return -1;

		}
	}

	

	////////////////////////////////////////////////////////
	//   Open file
	////////////////////////////////////////////////////////

	//Adding name extensions: file_name_integration_dt_cso_cn1m_cn2v_cn3t_stim_dur_stim_inc_min_max
	sprintf(file_ext,"%s_%.4f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f.asc",
		methods[integration].c_str(),dt,c_so,c_n1m,c_n2v,c_n3t,
		stim_dur,stim_inc,MIN_c,MAX_c);
	
	sprintf(file_spikes,"%s_spikes_%s",file_name,file_ext);
	sprintf(file_name,"%s_%s",file_name,file_ext);


	f = fopen(file_name,"w");
	f_spks = fopen(file_spikes,"w");

	if(!f|!f_spks)
	{
		cerr << "Error: error openning files"<<endl;
	}


	///////////////////////////////////////
	//Initializing Neurons and connections
	///////////////////////////////////////

	VavoulisModel so(VavoulisModel::SO);
	VavoulisModel n1m(VavoulisModel::N1M);
	VavoulisModel n2v(VavoulisModel::N2v);
	VavoulisModel n3t(VavoulisModel::N3t);

	std::vector<VavoulisModel *> neurons({&so,&n1m,&n2v,&n3t});
	std::vector<vector< VavoulisSynapse * > > syns(neurons.size());

	std::vector<double> prevs(neurons.size(),1); //Auxiliar vector for spike detection


	//Assigning each neuron a current value reference. 

	std::vector<double> c_values({c_so,c_n1m,c_n2v,c_n3t});
	
	//Auxiliar vector for non-feeding simulation
	std::vector<double> c_values_feed({c_so,0,c_n2v,25});
	std::vector<double> c_values_save;

	VavoulisSynapse n1m_so,n2v_so,so_n2v,n1m_n3t,n3t_n1m,n3t_n2v,n1m_n2v,n2v_n1m;

	//FAST 50 INHIB -90
	//SLOW 200 EXCIT 0
	if(connection >= 3) //All neurons connected
	{
		n1m_so.setSynapse(&n1m,&so,4.0,SLOW,EXCIT) ;
		syns[VavoulisModel::N1M].push_back(&n1m_so);

		n2v_so.setSynapse(&n2v,&so,1.0,SLOW,EXCIT);
		syns[VavoulisModel::N2v].push_back(&n2v_so);

		so_n2v.setSynapse(&so,&n2v,8.0,FAST,INHIB);
		syns[VavoulisModel::SO].push_back(&so_n2v);

	}

	if(connection >= 2) //N1M, N2v and N3t connected
	{
		n1m_n3t.setSynapse(&n1m,&n3t,8.0,FAST,INHIB) ;
		syns[VavoulisModel::N1M].push_back(&n1m_n3t);

		n3t_n1m.setSynapse(&n3t,&n1m,0.5,FAST,INHIB) ;
		syns[VavoulisModel::N3t].push_back(&n3t_n1m);

		n3t_n2v.setSynapse(&n3t,&n2v,2.0,FAST,INHIB) ;
		syns[VavoulisModel::N3t].push_back(&n3t_n2v);
	}
	if(connection >= 1) //N1M and N2v connected
	{
		n1m_n2v.setSynapse(&n1m,&n2v,50.0,FAST,INHIB) ;
		syns[VavoulisModel::N1M].push_back(&n1m_n2v);

		n2v_n1m.setSynapse(&n2v,&n1m,0.077,SLOW,EXCIT);
		syns[VavoulisModel::N2v].push_back(&n2v_n1m);

	}


	///////////////////////////////////////
	//Calculating iterations
	///////////////////////////////////////

	stim_dur*=1000; //To ms
	int stim_dur_iters = stim_dur/dt;

	if(secs_dur == -1) //If no time is specified, it is calculated based on the ramp paramaters.
	{
		//Compute number of iterations necessaries for the ramp
		secs_dur = ((MAX_c-MIN_c)/stim_inc)*rounds;
		iters = secs_dur*(stim_dur_iters);
	}
	else //If duration is specified, use it for the total iterations in the model. 
		iters = (secs_dur*1000 )/dt;

	feed_ini= (feed_ini*1000)/dt;
	feed_end= (feed_end*1000)/dt;
	cout << feed_ini << " " << feed_end<< endl;


	///////////////////////////////////////
	//Ramp initialization
	///////////////////////////////////////


	//Initialization of the Ramp Generator.
	RampGenerator rg(MIN_c,MAX_c,stim_inc,stim_dur);


	///////////////////////////////////////
	//Print parameters used. 
	///////////////////////////////////////

	printf("dt: %f \n",dt);
	printf("Ramp rounds: %d\n",rounds);
	printf("Simulation duration in ms: %.3f \n",secs_dur);
	printf("Iterations: %d \n",iters);
	printf("Parameters value\n");
	printf("File: %s\n",file_name);
	printf("File spikes: %s\n",file_spikes);
	printf("Connection %d\n",connection );
	printf("Integration method %s\n",methods[integration].c_str() );
	printf("Connection %d\n",connection );
	printf("-c_so=%.2f c_n1m=%.2f c_n2v=%.2f c_n3t=%.2f\n",c_so,c_n1m,c_n2v,c_n3t);
	printf("Ramp parameters:\n");
	printf("Min_c=%.2f Max_c=%.2f\n",MIN_c,MAX_c );
	printf("stim_dur=%.2f",stim_dur);
	printf(" stim_inc=%.2f\n",stim_inc);
	printf("feed_ini=%.2f feed_end=%.2f\n",feed_ini,feed_end );
	cout << endl;


	//Write File header 
	header = headers[connection].c_str();

	fprintf(f, "%s\n",header );
	fprintf(f_spks,"%f\n",SPIKE_TH);
	fprintf(f_spks, "%s\n",header );


	
	////////////////////////////////////////////////////////
	//   Simulation Mode Parameters
	////////////////////////////////////////////////////////

	int serie =0;//Variable used to reduce output file dimension. 
	double t=0.0;
	double c = MIN_c;

	//Starting clock
	clock_t begin = clock();

	for (int i=0; i < iters; i++)
	{
		//Write file each 3 iterations to reduce output file dimension
      serie = (serie + 1) % 4;
      if (serie == 3)
      {

		if(connection == 0)
			fprintf(f, "%f %f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V(),c);
		else if(connection == 1)
			fprintf(f, "%f %f %f\n", t,n1m.V(),n2v.V());
		else if(connection == 2)
			fprintf(f, "%f %f %f %f\n", t,n1m.V(),n2v.V(),n3t.V());
		else if(connection == 3)
			fprintf(f, "%f %f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V(),c);
		else if(connection == 4)
			fprintf(f, "%f %f %f %f %f %f %f %f %f\n", t,so.V(),so.getIsyn(),n1m.V(),n1m.getIsyn(),n2v.V(),n2v.getIsyn(),n3t.V(),n3t.getIsyn());
		else
			fprintf(f, "%f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V());

	}



		//////////////////////////////////////////////////////////////
		////////////////// SIMULATING SATIATED BEHAVIOUR /////////////
		//////////////////////////////////////////////////////////////

		if(feed_ini==i)
		{
			c_values_save = c_values; //save current values
			c_values = c_values_feed; //assign feeding current values
		}
		if(feed_end==i)
		{
			c_values = c_values_save; //restore current values
		}
		///////////////////////////////////////////////////////////

		//Update all parameters
		c = update_all(t,neurons,syns,integration,dt,c_values, &rg);

		//Detect spikes and write in spikes file.
  		detect_spikes(f_spks,neurons,prevs,t);

		t += dt;

		if(i==(int)iters/2)
			cout << "Half iterations performed" << endl;
	}


	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Execution time: %f\n",time_spent);

	printf("\n\n\n");

	fclose(f_spks);
	fclose(f);


	return 0;

}




double update_all(double _time, const vector<VavoulisModel *>& neurons,const std::vector<vector<VavoulisSynapse * > >& syns, integrators integr,double dt, const std::vector<double> & i_ext, RampGenerator * rg)
{
	if(integr == EULER)
	{
		update_euler(_time, dt, neurons, syns ,i_ext,rg);
	}
	else if(integr == RUNGE)
	{	
		update_runge(_time, neurons,syns, dt,i_ext, rg);	

	}

	for (int i=0; i< (int) neurons.size(); i++)
		if(i_ext[i]==-1)
			return rg->get_ext(i_ext[i],_time);


	return i_ext[VavoulisModel::N1M];

}


void update_runge(double _time,const std::vector<VavoulisModel * >& neurons,const std::vector<vector<VavoulisSynapse * > >& syns, double dt , const std::vector<double> &iext, RampGenerator * rg)
{
	intey(_time, neurons,syns, dt,iext, rg);

}



void funcion(double _time,const std::vector<std::vector<double > > & v_variables, std::vector<std::vector<double > >  &v_fvec, const std::vector<VavoulisModel *>& neurons,const std::vector<vector<VavoulisSynapse *> >& syns , const std::vector<double> & iext, RampGenerator * rg)
{	

	int n_vars=VavoulisModel::getNVars(); //Obtain number of variables in each neuron in the model
	int n_vars_syns=VavoulisSynapse::getNVars(); //Obtain number of of variables in each synapse in the model

	std::vector<double> ret(n_vars); //Auxiliar return vector for neuron variables
	std::vector<double> ret_syn(n_vars_syns); //Auxiliar return vector for synapse variables


	int pre_index =-1; //Presynaptic neuron type
	double vpre=-1; //Presynaptic neuron voltage value
	double i_syn=0; //Synaptic current 
	double i_ext =0; //Injected current 

	//Iterate through neurons array 
	for(int i=0; i<(int)neurons.size(); i++)
	{
		//reset return vector. 
		v_fvec[i].clear();
		i_syn=0;

		//Iterate through each neuron synapse
		for(int j=0; j<(int)syns[i].size(); j++)
		{
			int ref = j*n_vars_syns+n_vars;

			//Get update vector from the synapse 
			std::vector<double> vars_syns(v_variables[i].begin()+ref,v_variables[i].begin()+ref+n_vars_syns);

			//Obtain synaptic current
			i_syn += syns[i][j]->Isyn(vars_syns,v_variables[i][0]); //v_value
			
			//Get Vpre value as the Vs of the presynaptic neuron in the synapse. 
			pre_index = syns[i][j]->getPreType();
			vpre = v_variables[pre_index][0];

			syns[i][j]->funcion_simple(_time, vars_syns, ret_syn, vpre);
			//Add syn variables to return vector. 
			v_fvec[i].insert(v_fvec[i].end(),ret_syn.begin(),ret_syn.end());
		}

		//Get update vector from the synapse. 
		std::vector<double> vars_neu(v_variables[i].begin(),v_variables[i].begin()+n_vars);
		i_ext = rg->get_ext(iext[i],_time);
		neurons[i]->funcion_simple(_time, vars_neu,ret, i_ext,i_syn);

		//Add neuron variables to return vector. 
		v_fvec[i].insert(v_fvec[i].begin(),ret.begin(),ret.end());

	}

}


/*======================================*/
/* Rutina de integración                */
/*======================================*/
double intey(double _time, const std::vector<VavoulisModel *> &neurons,const std::vector<vector<VavoulisSynapse *> >& syns, double inc_integracion, const std::vector<double> & iext, RampGenerator * rg)
{	
	std::vector<std::vector<double > > v_apoyo(neurons.size());
	std::vector<std::vector<double > > v_variables(neurons.size());
	std::vector<std::vector<double > > v_retorno(neurons.size());
	double v_k[neurons.size()][6][MAX_VARS];
	int n_total_variables[N_NEU];

	double u=0.0;
	int j;
	int i;
	int total_variables;

	for(i=0; i < (int)neurons.size(); ++i)
	{
		total_variables = neurons[i]->getNVars()+syns[i].size()*VavoulisSynapse::getNVars();
		n_total_variables[i]=total_variables;

		std::vector<double> vars_neu = neurons[i]->getVariables();
		v_variables[i].insert(v_variables[i].end(),vars_neu.begin(),vars_neu.end());

		for(int j=0; j<(int)syns[i].size(); j++)
		{
			//Get update vector from the synapse 
			std::vector<double> vars_syns = syns[i][j]->getVariables();
			v_variables[i].insert(v_variables[i].end(),vars_syns.begin(),vars_syns.end());
		}

	}

	funcion(_time, v_variables,v_retorno,neurons,syns,iext, rg);


	for(i=0; i < (int)neurons.size(); ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
		{
		  	v_k[i][0][j]=inc_integracion*v_retorno[i][j];
		  	//Update neurons
		  	v_apoyo[i].resize(v_variables[i].size());
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*.2;
		} 
	}




	funcion(_time+inc_integracion/5, v_apoyo,v_retorno,neurons,syns,iext, rg);

	for(i=0; i < (int)neurons.size(); ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
		{
		  	v_k[i][1][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*.075+v_k[i][1][j]*0.225;
		} 
	}

	funcion(_time+inc_integracion*0.3, v_apoyo,v_retorno,neurons,syns,iext, rg);

	for(i=0; i < (int)neurons.size(); ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
	  	{
		  	v_k[i][2][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*.3-v_k[i][1][j]*0.9+v_k[i][2][j]*1.2;
		}
	}

	funcion(_time+inc_integracion*0.6, v_apoyo,v_retorno,neurons,syns,iext,  rg);

	for(i=0; i < (int)neurons.size(); ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
	  	{
		  	v_k[i][3][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*0.075+v_k[i][1][j]*0.675-v_k[i][2][j]*0.6+v_k[i][3][j]*0.75;
	  	} 	
	}

	funcion(_time+inc_integracion*0.9, v_apoyo,v_retorno,neurons,syns,iext,  rg);
	for(i=0; i < (int)neurons.size(); ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
		{
		  	v_k[i][4][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*0.660493827160493
		       +v_k[i][1][j]*2.5
		       -v_k[i][2][j]*5.185185185185185
		       +v_k[i][3][j]*3.888888888888889
		       -v_k[i][4][j]*0.864197530864197;
		}
	}

	funcion(_time+inc_integracion, v_apoyo,v_retorno,neurons,syns,iext, rg);

	for(i=0; i < (int)neurons.size(); ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
		{
		  	v_k[i][5][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*0.1049382716049382+
			       v_k[i][2][j]*0.3703703703703703+
			       v_k[i][3][j]*0.2777777777777777+
			       v_k[i][4][j]*0.2469135802469135;
		}

	 }
		
	for(i=0; i < (int)neurons.size(); ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
		{
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*0.098765432098765+
	                       v_k[i][2][j]*0.396825396825396+
	                       v_k[i][3][j]*0.231481481481481+
	                       v_k[i][4][j]*0.308641975308641-
	                       v_k[i][5][j]*0.035714285714285;
		}

	 }
	
	for(i=0; i < (int)neurons.size(); ++i)
	{
		neurons[i]->set_variables(v_apoyo[i]);
		int n= neurons[i]->getNVars();

	  	for(int j=0; j<(int)syns[i].size(); j++)
		{
			//Get update vector from the synapse 
			std::vector<double> vars_syns(v_apoyo[i].begin()+j*2+n,v_apoyo[i].begin()+j*2+n+2);
			syns[i][j]->set_variables(vars_syns);
		}
	 }
	
  return u;
}


void update_euler(double _time, double dt,std::vector<VavoulisModel *> neurons, std::vector<std::vector<VavoulisSynapse* > > syns ,const std::vector<double> &iext,RampGenerator * rg)
{
	double isyn =0;
	double i_ext=0;

	for(int i=0; i<(int)neurons.size(); i++)
	{
		isyn=0;
		for(int j=0; j< (int)syns[i].size();j++)
		{
			isyn+=syns[i][j]->Isyn(); //Obtain the synaptic current
			syns[i][j]->update_variables(dt,_time); //Update synapse variables
		}
		i_ext = rg->get_ext(iext[i],_time); 
		neurons[i]->update_variables( dt, _time, i_ext,isyn);
	}

}



void detect_spikes(FILE * f_spks, std::vector<VavoulisModel *> const & neurons, std::vector<double> &prevs, double t )
{

	bool fst_inrow = true;
	string buff ="";
	double dev;

	for(int n=0; n<(int)neurons.size(); ++n)
	{	
		dev = neurons[n]->getdV();
		//0.001 less than that change in the derivative is not a spike
		//NOTE: it might be necessary to adjust MIN_SPIKE_CHANGE for different spike shapes. 
		if(prevs[n] >0 && dev <0 && (prevs[n]-dev)>MIN_SPIKE_CHANGE) //If it's a spike (from pos derivate to 0 derivate)
		{
			fst_inrow = false;

			if(neurons[n]->V() >SPIKE_TH) //Spike must be over a threshold. 
				buff += to_string(neurons[n]->V()) +" ";
		}
		else
				buff +=", "; //No spike value. 

		prevs[n]=dev;

	}
	if(!fst_inrow)
	{
			fprintf(f_spks,"%f %s\n",t,buff.c_str());
	}
		
}


int parse_input(int argc,char *argv[], void ** arguments)
{

	int * aux_i;
	char ** aux_s;
	double * aux_d;
	float * aux_f;

	int num_args = *(&arg_names + 1) - arg_names;

	for(int i=1; i<argc; i+=2)
	{
		//Look for the argument
		for(int j=0; j<num_args; j++)
		{
			if(argv[i][0]!='-')
			{
				break;
			}

			if( strcmp(argv[i], arg_names[j].c_str()) == 0)
			{
				switch(arg_types[j])
				{
					case Integer:
						aux_i = (int *) arguments[j];
						*aux_i = atoi(argv[i+1]);
						break;
	
					case String:
						aux_s =(char **) arguments[j];
						*aux_s =argv[i+1];
						break;
	
					case IntegrationMeth:
						aux_i = (int *) arguments[j];
						if (strcmp(argv[i+1], "-e") == 0){
							*aux_i = 0;
						}
						else if(strcmp(argv[i+1], "-r") == 0)
							*aux_i = 1;

						break;
	
					case Float:
						aux_f = (float *) arguments[j];
						*aux_f = atof(argv[i+1]);
						break;
					
					case Double:
						aux_d = (double *) arguments[j];
						*aux_d = atof(argv[i+1]);
						break;
					
					default:
						cout << "Argument type not defined" << endl;
						return ERROR;

				}
				break;
			}
			else if(j == num_args-1)
			{
				cerr << "Incorrect argument key: " << argv[i]  << endl;
						return ERROR;
			}

		}
	}

	return OK;

}

void show_help()
{

	cout <<format << endl;
	
	cout << "-connection: specifies connection type in the circuit" << endl;
	cout << "\t 0 all neurons isolated" << endl;
	cout << "\t 1 N1M and N2v are connected" << endl;
	cout << "\t 2 N1M, N2v and N3t are connected" << endl;
	cout << "\t 3 complete circuit: N1M, N2v, N3t and SO are connected" << endl;
	cout << endl;
	cout << "-file_name: name of the file where V info will be written. Spikes will be recorded in file_name_spikes"<< endl;
	cout << endl;
	cout << "-integration_method:"<<endl;
	cout << "-e for Euler"<<endl;
	cout << "-r for Runge-Kutta"<<endl;
	cout << endl;
	cout << "-c_so/c_n1m/c_n2v/c_n3t: current values applied to each neuron respectivelly"<< endl;
	cout << "default values: 10 6 4 0"<<endl;
	cout << "when value is -1 the applied current is the ramp generated" << endl;
	cout << endl;
	cout << "secs_dur: Force specific duration of the simulation"<<endl;
	cout << "\t IMPORTANT: If this parameter is omitted the duration will be computed for 2 up-down ramps performance"<<endl;
	cout << "Ramp parameters:"<<endl;
	cout << "\t stim_dur: Time that the ramp stays in the same value"<<endl;
	cout << "\t stim_inc: Current increment each stim_dur"<< endl;
	cout << "\t MIN_c: minimum current value"<<endl;
	cout << "\t MAX_c: maximum current value"<<endl;
	cout << endl;
	cout << "\t rounds: number of rounds in the ramp."<<endl;
	cout << "\t\t one round is going from min to max."<<endl;
	cout << endl;

}
