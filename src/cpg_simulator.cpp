/*************************************************************
	Developed by Alicia Garrido Peña (2020)

	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	
	Please, if you use this implementation cite the two papers above in your work. 
*************************************************************/


#include "cpg_simulator.h"


#include <iostream>
using namespace std; 

CPGSimulator::CPGSimulator()
{
	n_neurons=0;
	connection=-1;
}

CPGSimulator::CPGSimulator(int connection, std::vector<double> c_values, RampGenerator rg)
{
	init(connection,c_values,rg);
}

int CPGSimulator::init(int connection, std::vector<double> c_values, RampGenerator rg)
{


	///////////////////////////////////////
	//Initializing Neurons and connections
	///////////////////////////////////////

	if(connection < 0)return 0;

	this->connection=connection;

	VavoulisModel so(VavoulisModel::SO);
	VavoulisModel n1m(VavoulisModel::N1M);
	VavoulisModel n2v(VavoulisModel::N2v);
	VavoulisModel n3t(VavoulisModel::N3t);

	neurons.assign({so,n1m,n2v,n3t});
	n_neurons=neurons.size();
	syns.resize(n_neurons);

	//Assigning each neuron a current value reference. 

	this->c_values = c_values;

	VavoulisSynapse n1m_so,n2v_so,so_n2v,n1m_n3t,n3t_n1m,n3t_n2v,n1m_n2v,n2v_n1m;


	//FAST 50 INHIB -90
	//SLOW 200 EXCIT 0
	if(connection >= 3) //All neurons connected
	{
		n1m_so.setSynapse(&neurons[VavoulisModel::N1M],&neurons[VavoulisModel::SO],4.0,SLOW,EXCIT) ;
		syns[VavoulisModel::N1M].push_back(n1m_so);

		n2v_so.setSynapse(&neurons[VavoulisModel::N2v],&neurons[VavoulisModel::SO],1.0,SLOW,EXCIT);
		syns[VavoulisModel::N2v].push_back(n2v_so);

		so_n2v.setSynapse(&neurons[VavoulisModel::SO],&neurons[VavoulisModel::N2v],8.0,FAST,INHIB);
		syns[VavoulisModel::SO].push_back(so_n2v);

	}

	if(connection >= 2) //N1M, N2v and N3t connected
	{
		n1m_n3t.setSynapse(&neurons[VavoulisModel::N1M],&neurons[VavoulisModel::N3t],8.0,FAST,INHIB) ;
		syns[VavoulisModel::N1M].push_back(n1m_n3t);

		n3t_n1m.setSynapse(&neurons[VavoulisModel::N3t],&neurons[VavoulisModel::N1M],0.5,FAST,INHIB) ;
		syns[VavoulisModel::N3t].push_back(n3t_n1m);

		n3t_n2v.setSynapse(&neurons[VavoulisModel::N3t],&neurons[VavoulisModel::N2v],2.0,FAST,INHIB) ;
		syns[VavoulisModel::N3t].push_back(n3t_n2v);
	}
	if(connection >= 1) //N1M and N2v connected
	{
		n1m_n2v.setSynapse(&neurons[VavoulisModel::N1M],&neurons[VavoulisModel::N2v],50.0,FAST,INHIB) ;
		syns[VavoulisModel::N1M].push_back(n1m_n2v);

		n2v_n1m.setSynapse(&neurons[VavoulisModel::N2v],&neurons[VavoulisModel::N1M],0.077,SLOW,EXCIT);
		syns[VavoulisModel::N2v].push_back(n2v_n1m);
	}


	///////////////////////////////////////
	//Ramp initialization
	///////////////////////////////////////


	this->rg = rg;

	return 1;

}

void CPGSimulator::print()
{
	for (int i=0;i<n_neurons;i++)
	{
		cout << "Neuron: ";
		neurons[i].print();
		cout << "\nSynapses: \n\t";
		cout << "Number of synapses: "<< syns[i].size()<< endl;
		for(int j=0;j<(int)syns[i].size();j++)
		{
			cout<< "\tSynapse number: "<< j<<endl;
			syns[i][j].print();
			cout << endl;
		}
	}
}


void CPGSimulator::write(FILE *f,double t,double c)
{
	if(connection == 0 || connection==3)
		fprintf(f, "%f %f %f %f %f %f\n",
		 t,neurons[VavoulisModel::SO].V(),neurons[VavoulisModel::N1M].V(),neurons[VavoulisModel::N2v].V(),neurons[VavoulisModel::N3t].V(),c);
	else if(connection == 1)
		fprintf(f, "%f %f %f\n", t,neurons[VavoulisModel::N1M].V(),neurons[VavoulisModel::N2v].V());
	else if(connection == 2)
		fprintf(f, "%f %f %f %f\n", t,neurons[VavoulisModel::N1M].V(),neurons[VavoulisModel::N2v].V(),neurons[VavoulisModel::N3t].V());
	else if(connection == 4)
		fprintf(f, "%f %f %f %f %f %f %f %f %f\n", 
		t,neurons[VavoulisModel::SO].V(),neurons[VavoulisModel::SO].getIsyn(),neurons[VavoulisModel::N1M].V(),neurons[VavoulisModel::N1M].getIsyn(),neurons[VavoulisModel::N2v].V(),neurons[VavoulisModel::N2v].getIsyn(),neurons[VavoulisModel::N3t].V(),neurons[VavoulisModel::N3t].getIsyn());
	else
		fprintf(f, "%f %f %f %f %f\n", t,neurons[VavoulisModel::SO].V(),neurons[VavoulisModel::N1M].V(),neurons[VavoulisModel::N2v].V(),neurons[VavoulisModel::N3t].V());

}


void CPGSimulator::simulate(FILE * f,FILE * f_spks,double iters,double dt,integrators integration,double satiated_ini,double satiated_end)
{
	
	////////////////////////////////////////////////////////
	//   Simulation Mode Parameters
	////////////////////////////////////////////////////////

	
	int serie =0;//Variable used to reduce output file dimension. 
	double t=0.0;
	double c =0;

	std::vector<double> prevs(n_neurons,1);//Auxiliar vector for spike detection

	// "In satiated animals, N3t keeps the feeding network under its suppressive control." [1]
	// To simulate satiated activity N3t is stimulated by the injected current, while N1M value is supressed.
	std::vector<double> c_values_staited({c_values[VavoulisModel::SO],0,c_values[VavoulisModel::N2v],25});
	std::vector<double> c_values_save;

	for (int i=0; i < iters; i++)
	{

      serie = (serie + 1) % 4;
      if (serie == 3)
      {
		write(f,t,c);

	  }

		//////////////////////////////////////////////////////////////
		/////////////// SIMULATING SATIATED BEHAVIOUR ////////////////
		//////////////////////////////////////////////////////////////

		if(satiated_ini==i)
		{
			c_values_save = c_values; //save current values
			c_values = c_values_staited; //assign satiated current values
		}
		if(satiated_end==i)
		{
			c_values = c_values_save; //restore current values
		}


		///////////////////////////////////////////////////////////

		//Integrate variables in the model. 
		c = update_all(t,integration,dt);

		//Detect spikes and write in spikes file.
  		detect_spikes(f_spks,prevs,t);

		t += dt;

		if(i==(int)iters/2)
			cout << "Half iterations" << endl;
	}


}


void CPGSimulator::detect_spikes(FILE * f_spks, std::vector<double> &prevs, double t )
{

	bool fst_inrow = true;
	string buff ="";
	double dev;

	for(int n=0; n<n_neurons; ++n)
	{	
		dev = neurons[n].getdV();
		//0.001 less than that change in the derivative is not a spike
		//NOTE: it might be necessary to adjust MIN_SPIKE_CHANGE for different spike shapes. 
		if(prevs[n] >0 && dev <0 && (prevs[n]-dev)>MIN_SPIKE_CHANGE) //If it's a spike (from pos derivate to 0 derivate)
		{
			fst_inrow = false;

			if(neurons[n].V() >SPIKE_TH)
				buff += to_string(neurons[n].V()) +" ";
		}
		else
				buff +=", ";

		prevs[n]=dev;

	}
	if(!fst_inrow)
	{
			fprintf(f_spks,"%f %s\n",t,buff.c_str());
	}
		
}



double CPGSimulator::update_all(double _time, integrators integr,double dt)
{
	if(integr == EULER)
	{
		update_euler(_time, dt);
	}
	else if(integr == RUNGE)
	{	
		update_runge(_time, dt);	

	}


	for (int i=0; i< n_neurons; i++)
		if(c_values[i]==-1)
			return rg.get_ext(c_values[i],_time);


	return c_values[VavoulisModel::N1M];

}


void CPGSimulator::update_euler(double _time, double dt)
{
	double isyn =0;
	double i_ext=0;

	for(int i=0; i<n_neurons; i++)
	{
		isyn=0;
		for(int j=0; j< (int)syns[i].size();j++)
		{
			isyn+=syns[i][j].Isyn();
			syns[i][j].update_variables(dt,_time);
		}
		i_ext = rg.get_ext(c_values[i],_time);
		neurons[i].update_variables( dt, _time, i_ext,isyn);
	}

}


void CPGSimulator::update_runge(double _time, double dt)
{
	intey(_time, dt);

}


void CPGSimulator::diffs(double _time,const std::vector<std::vector<double > > & v_variables, std::vector<std::vector<double > >  &v_fvec)
{	

	int n_vars=VavoulisModel::getNVars();
	int n_vars_syns=VavoulisSynapse::getNVars();

	std::vector<double> ret(n_vars);
	std::vector<double> ret_syn(n_vars_syns);


	int pre_index =-1;
	double vpre=-1;
	double i_syn=0;
	double i_ext =0;

	//Iterate through neurons array 
	for(int i=0; i<n_neurons; i++)
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
			i_syn += syns[i][j].Isyn(vars_syns,v_variables[i][0]); //v_value
			
			//Get Vpre value as the Vs of the presynaptic neuron in the synapse. 
			pre_index = syns[i][j].getPreType();
			vpre = v_variables[pre_index][0];

			syns[i][j].diffs_fun(_time, vars_syns, ret_syn, vpre);
			//Add syn variables to return vector. 
			v_fvec[i].insert(v_fvec[i].end(),ret_syn.begin(),ret_syn.end());
		}

		//Get update vector from the synapse. 
		std::vector<double> vars_neu(v_variables[i].begin(),v_variables[i].begin()+n_vars);
		i_ext = rg.get_ext(c_values[i],_time);
		neurons[i].diffs_fun(_time, vars_neu,ret, i_ext,i_syn);

		//Add neuron variables to return vector. 
		v_fvec[i].insert(v_fvec[i].begin(),ret.begin(),ret.end());

	}

}


/*======================================*/
/* Rutina de integración                */
/*======================================*/
double CPGSimulator::intey(double _time, double inc_integracion)
{	
	std::vector<std::vector<double > > v_apoyo(n_neurons);
	std::vector<std::vector<double > > v_variables(n_neurons);
	std::vector<std::vector<double > > v_retorno(n_neurons);
	double v_k[n_neurons][6][MAX_VARS];
	int n_total_variables[n_neurons];

	double u=0.0;
	int j;
	int i;
	int total_variables;

	for(i=0; i < n_neurons; ++i)
	{
		total_variables = neurons[i].getNVars()+syns[i].size()*VavoulisSynapse::getNVars();
		n_total_variables[i]=total_variables;

		std::vector<double> vars_neu = neurons[i].getVariables();
		v_variables[i].insert(v_variables[i].end(),vars_neu.begin(),vars_neu.end());

		for(int j=0; j<(int)syns[i].size(); j++)
		{
			//Get update vector from the synapse 
			std::vector<double> vars_syns = syns[i][j].getVariables();
			v_variables[i].insert(v_variables[i].end(),vars_syns.begin(),vars_syns.end());
		}

	}

	diffs(_time, v_variables,v_retorno);


	for(i=0; i < n_neurons; ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
		{
		  	v_k[i][0][j]=inc_integracion*v_retorno[i][j];
		  	//Update neurons
		  	v_apoyo[i].resize(v_variables[i].size());
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*.2;
		} 
	}




	diffs(_time+inc_integracion/5, v_apoyo,v_retorno);

	for(i=0; i < n_neurons; ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
		{
		  	v_k[i][1][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*.075+v_k[i][1][j]*0.225;
		} 
	}

	diffs(_time+inc_integracion*0.3, v_apoyo,v_retorno);

	for(i=0; i < n_neurons; ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
	  	{
		  	v_k[i][2][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*.3-v_k[i][1][j]*0.9+v_k[i][2][j]*1.2;
		}
	}

	diffs(_time+inc_integracion*0.6, v_apoyo,v_retorno);

	for(i=0; i < n_neurons; ++i)
	{
		for(j=0;j<n_total_variables[i];++j)
	  	{
		  	v_k[i][3][j]=inc_integracion*v_retorno[i][j];
			v_apoyo[i][j]=v_variables[i][j]+v_k[i][0][j]*0.075+v_k[i][1][j]*0.675-v_k[i][2][j]*0.6+v_k[i][3][j]*0.75;
	  	} 	
	}

	diffs(_time+inc_integracion*0.9, v_apoyo,v_retorno);
	for(i=0; i < n_neurons; ++i)
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

	diffs(_time+inc_integracion, v_apoyo,v_retorno);

	for(i=0; i < n_neurons; ++i)
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
		
	for(i=0; i < n_neurons; ++i)
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
	
	for(i=0; i < n_neurons; ++i)
	{
		neurons[i].set_variables(v_apoyo[i]);
		int n= neurons[i].getNVars();

	  	for(int j=0; j<(int)syns[i].size(); j++)
		{
			//Get update vector from the synapse 
			std::vector<double> vars_syns(v_apoyo[i].begin()+j*2+n,v_apoyo[i].begin()+j*2+n+2);
			syns[i][j].set_variables(vars_syns);
		}
	 }
	
  return u;
}


