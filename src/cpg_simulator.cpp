int CPGSimulator::init()
{


	///////////////////////////////////////
	//Initializing Neurons and connections
	///////////////////////////////////////


	VavoulisModel so(VavoulisModel::SO);
	VavoulisModel n1m(VavoulisModel::N1M);
	VavoulisModel n2v(VavoulisModel::N2v);
	VavoulisModel n3t(VavoulisModel::N3t);

	// std::vector<VavoulisModel *> neurons({&so,&n1m,&n2v,&n3t});
	std::vector<VavoulisModel> neurons({so,n1m,n2v,n3t});

	std::vector<vector< VavoulisSynapse * > > syns(neurons.size());



	std::vector<double> prevs(neurons.size(),1);


	//Assigning each neuron a current value reference. 

	this->c_values = c_values;
	//"When feeding a high current is applied to N3t and low to N1M so it ceases its activity"
	//it is suppresed in N3t to achieve tonic spiking
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
	//Ramp initialization
	///////////////////////////////////////


	//Initialization of the Ramp Generator.
	RampGenerator rg(MIN_c,MAX_c,stim_inc,stim_dur);
	this->rg = rg;


}

CPGSimulator::simulate()
{
	
	////////////////////////////////////////////////////////
	//   Simulation Mode Parameters
	////////////////////////////////////////////////////////

	//Variable used to reduce output file dimension. 
	int serie =0;
	double t=0.0;
	double c = MIN_c;

	for (int i=0; i < iters; i++)
	{

 //      serie = (serie + 1) % 4;
 //      if (serie == 3)
 //      {

		// fprintf(f, "%f %f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V(),c);
		if(connection == 0)
			fprintf(f, "%f %f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V(),c);
		else if(connection == 1)
			// fprintf(f, "%f %f %f\n", t,n1m.V(),n2v.V());
			fprintf(f, "%f %f %f %f %f\n", t,n1m.V(),n2v.V(),n1m.getIsyn(),n2v.getIsyn());
		else if(connection == 2)
			fprintf(f, "%f %f %f %f\n", t,n1m.V(),n2v.V(),n3t.V());
		else if(connection == 3)
			fprintf(f, "%f %f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V(),c);
		else if(connection == 4)
			// fprintf(f, "%f %f %f %f %f %f %f %f %f %f\n", t,so.V(),so.publics[VavoulisModel::_isyn,n1m.V(),n1m.publics[VavoulisModel::_isyn,n2v.V(),n2v.publics[VavoulisModel::_isyn,n3t.V(),n3t.publics[VavoulisModel::_isyn,c);
			// fprintf(f, "%f %f %f %f %f %f %f %f %f\n", t,so.V(),so.Va(),n1m.V(),n1m.Va(),n2v.V(),n2v.Va(),n3t.V(),n3t.Va());
		fprintf(f, "%f %f %f %f %f %f %f %f %f\n", t,so.V(),so.getIsyn(),n1m.V(),n1m.getIsyn(),n2v.V(),n2v.getIsyn(),n3t.V(),n3t.getIsyn());
		else
			fprintf(f, "%f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V());

	// }



		//////////////////////////////////////////////////////////////
		////////////////// SIMULATING FEEDING ////////////////////////
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



		// c = current;
		//10,6,2,0 
		//			<- con 4.4 en n2v no va
		//simulating feeding 10 8 2 0 con 15 dentro del bucle
		//N1-driven 8.5,c,2,0
		//SO-driven c,10,1,4 <- sin 4 en n3 no varía.

		//SO no curr 0 6 2 0

		c = update_all(t,neurons,syns,integration,dt,c_values, &rg);

  		detect_spikes(f_spks,neurons,prevs,t);

		t += dt;

		if(i==(int)iters/2)
			cout << i << endl;
	}


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
	// update_synapses(integr,dt);
	// return -1;

}


double update_runge(double _time,const std::vector<VavoulisModel * >& neurons,const std::vector<vector<VavoulisSynapse * > >& syns, double dt , const std::vector<double> &iext, RampGenerator * rg)
{
	// intey(_time, neurons, dt,iext, rg);
	intey(_time, neurons,syns, dt,iext, rg);
	// return iext;
	return 0;

}



void funcion(double _time,const std::vector<std::vector<double > > & v_variables, std::vector<std::vector<double > >  &v_fvec, const std::vector<VavoulisModel *>& neurons,const std::vector<vector<VavoulisSynapse *> >& syns , const std::vector<double> & iext, RampGenerator * rg)
{	

	int n_vars=VavoulisModel::getNVars();
	int n_vars_syns=VavoulisSynapse::getNVars();

	std::vector<double> ret(n_vars);
	std::vector<double> ret_syn(n_vars_syns);


	int pre_index =-1,n;
	double vpre=-1;
	double i_syn=0;
	double i_ext =0;

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


double update_euler(double _time, double dt,std::vector<VavoulisModel *> neurons, std::vector<std::vector<VavoulisSynapse* > > syns ,const std::vector<double> &iext,RampGenerator * rg)
{
	double isyn =0;
	double i_ext=0;

	for(int i=0; i<(int)neurons.size(); i++)
	{
		isyn=0;
		for(int j=0; j< (int)syns[i].size();j++)
		{
			isyn+=syns[i][j]->Isyn();
			syns[i][j]->update_variables(dt,_time);
		}
		i_ext = rg->get_ext(iext[i],_time);
		neurons[i]->update_variables( dt, _time, i_ext,isyn);
	}

	return iext[0];
}

