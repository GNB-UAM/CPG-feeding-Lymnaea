/*
	Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
	and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
	Please, if you use this implementation cite the two paper above in your work. 
*/

/*
./lymnaea-va -connection 3 -file_name ./data/val -integrator -r -c_so 10 -c_n1m 6 -c_n2v 2 -c_n3t 0 -dt 0.001 -stim_dur 2 -stim_inc 5 -MIN_c 0 -MAX_c 10 -secs_dur 4

make vavocplus connection=3 file_name=prueba_time integration=-r dt=0.0006 c_so=10 c_n1m=6 c_n2v=2 c_n3t=0 stim_dur=4.6 stim_inc=0.5 MIN_c=0 MAX_c=10 iters=9.2

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <string>
#include <vector>
#include <iostream>

// #include "vavoulis_neuron.h"
// #include "ramp_generator.h"
// #include "vavoulis_synapse.h"
#include "cpg_simulator.h"

using namespace std;

#define n_variables 8

#define DUR_BURST 2300 //mseconds	
// #define MAX_VARS 18

// #define N_NEU 4
// #define N_VARS 6


// #define SPIKE_TH -50.0
// #define MIN_SPIKE_CHANGE 0.001



#define MAX_STRING 20000

// #define NUM_ARG 16
#define ERROR 0
#define OK 1

enum prim_types{String,Integer, Float, Double, IntegrationMeth};
string arg_names[]={"-connection","-file_name","-integrator","-dt","-c_so","-c_n1m","-c_n2v","-c_n3t","-stim_dur","-stim_inc","-MIN_c","-MAX_c","-secs_dur","-rounds","-feed_ini","-feed_end"};
string format = "Format: ./lymn -connection val -file_name val -integration_method -dt val val -c_so val -c_n1m val -c_n2v val -c_n3t val -stim_dur val -stim_inc val -MIN_c val -MAX_c val [-secs_dur val] -rounds val -feed_ini val -feed_end val\n";
prim_types arg_types[]={Integer,String,IntegrationMeth,Double,Double,Double,Double,Double,Double,Double,Double,Double,Double,Integer,Double,Double};

string headers[] = {"t SO N1M N2v N3t c", "t N1M N2v IsynN1M IsynN2v","t N1M N2v N3t","t SO N1M N2v N3t c","t SO IsynSO N1M IsynN1M N2v IsynN2v N3t IsynN3t",
			"t SO N1M N2v N3t", "t SO N1M N2v N3t c"};

//TODO: methods in new class?
string methods[] = {"Euler","Runge-Kutta"}; 


int parse_input(int argc,char *argv[], void ** arguments);
void show_help();

int main(int argc, char * argv[])
{

	int connection =3;
	char * file_name = "./data/complete.asc";
	char file_spikes[MAX_STRING];
	char file_ext[MAX_STRING];
	double secs_dur = -1;
	int iters =-1;
	CPGSimulator::integrators integration = CPGSimulator::EULER;
	double dt=0.01;
	double c_so=-1,c_n1m=-1,c_n2v=-1,c_n3t=-1;
	double stim_dur=-1;
	double stim_inc=-1,MIN_c=-1,MAX_c=-1;
	double feed_ini =  -1;
	double feed_end =  -1;

	int rounds=4;

	FILE *f,*f_spks;

	const char * header;

	///////////////////////////////////////
	//Input Parameters 
	///////////////////////////////////////

	if(argc == 1){
		cout <<format<< endl;
		cout << "If any of those parameters is skiped, the default value will be assigned\n ./lymn --help for more info" << endl;
		// printf("%d\n", argc);
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
		if(parse_input(argc,argv,arguments)==ERROR)return -1;
		// cout <<format<< endl;
	}

	

	////////////////////////////////////////////////////////
	//   Open file
	////////////////////////////////////////////////////////

	//Add currents
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
	//Calculating iterations
	///////////////////////////////////////

	stim_dur*=1000; //To ms
	int stim_dur_iters = stim_dur/dt;

	if(secs_dur == -1)
	{
		//Compute number of iterations necessaries for the ramp
		secs_dur = ((MAX_c-MIN_c)/stim_inc)*rounds;
		iters = secs_dur*(stim_dur_iters);
	}
	else //If duration is specified, use it for iters computation. 
		iters = (secs_dur*1000 )/dt;

	if(feed_ini >0 and feed_end >0){

		feed_ini= (feed_ini*1000)/dt;
		feed_end= (feed_end*1000)/dt;
	}
	// cout << feed_ini << " " << feed_end<< endl;


	///////////////////////////////////////
	//Ramp initialization
	///////////////////////////////////////


	//Initialization of the Ramp Generator.
	RampGenerator rg(MIN_c,MAX_c,stim_inc,stim_dur);

	//Write File header 
	header = headers[connection].c_str();

	fprintf(f, "%s\n",header );
	fprintf(f_spks,"%f\n",SPIKE_TH);
	fprintf(f_spks, "%s\n",header );


	CPGSimulator cpg;

	std::vector<double> c_values({c_so,c_n1m,c_n2v,c_n3t});

	cpg.init(connection,c_values);

	cpg.print();

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
	printf(" stim_inc=%.2f",stim_inc);
	cout << endl;




	//Starting clock
	clock_t begin = clock();


	cpg.simulate(f,f_spks,iters,dt,integration,feed_ini,feed_end,connection);
	

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Execution time: %f\n",time_spent);

	printf("\n\n\n");

	fclose(f_spks);
	fclose(f);


	return 0;

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


// int write_file(FILE * f,int connection)
// {
// 	if(connection == 0)
// 		fprintf(f, "%f %f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V(),c);
// 	else if(connection == 1)
// 		// fprintf(f, "%f %f %f\n", t,n1m.V(),n2v.V());
// 		fprintf(f, "%f %f %f %f %f\n", t,n1m.V(),n2v.V(),n1m_isyn,n2v_isyn);
// 	else if(connection == 2)
// 		fprintf(f, "%f %f %f %f\n", t,n1m.V(),n2v.V(),n3t.V());
// 	else if(connection == 3)
// 		fprintf(f, "%f %f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V(),c);
// 	// else if(connection == 4)
// 	// 	// fprintf(f, "%f %f %f %f %f %f %f %f %f %f\n", t,so.V(),so.publics[VavoulisModel::_isyn,n1m.V(),n1m.publics[VavoulisModel::_isyn,n2v.V(),n2v.publics[VavoulisModel::_isyn,n3t.V(),n3t.publics[VavoulisModel::_isyn,c);
// 		// fprintf(f, "%f %f %f %f %f %f %f %f %f\n", t,so.V(),so_isyn,n1m.V(),n1m_isyn,n2v.V(),n2v_isyn,n3t.V(),n3t_isyn);
// 	else
// 		fprintf(f, "%f %f %f %f %f\n", t,so.V(),n1m.V(),n2v.V(),n3t.V());
// }