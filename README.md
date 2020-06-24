# Feeding CPG Model Lymnaea Stagnalis #

This model is an implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x [1] and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020).	Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020. [2]

Please, if you use this implementation cite the two papers above in your work. 

You can find a summary of the model here implemented in [2] and a detailed description of the model and all the equations here implemented in [1].

A Central Pattern Generator (CPG) is a neural circuit which is capable to generate authonomously, robust sequencies of neural activity which are usually related to motor activity in many individuals, including humans.The model here implemented represents a feeding CPG in the mollusk Lymnaea Stagnalis, composed by its three main neurons leading the rhythm: N1M, N2v and N3t and a fourth modulator neuron called SO.

This CPG model is based in Hodking-Huxley formalism and is formed by two compartments model (soma-axon) with the slow and fast activity differentiated en each of them. The connections in the circuits are implemented by a gradual synapse also defined by Vavoulis et al. where the wave form in the activity has a key role. 

The neurons in the model are not endogenously variable, but it is possible to induce variability by changing the current injected to each neuron throghout the simulation. For this purpose, RampGenerator calculates stimulation ramp as used in [2], which can be applied to any of the four neurons in the circuit. 


## Configuration and first run

### Generate documentation
By running 
	
	make doxygen 

the documentation for this implementation will be generated in html and latex form. 
	
### Build 
If it is the first time you use this model first run:
	
	make init

This will create data and images directories. 

To build the project do:
	
	make 

This will create the executable file called feeding_cpg. By running 

	./feeding_cpg --help 

an explanation of the input arguments will prompt.

### Data recording

### Choosing integrator and integration increment
-integrator flag is used to choose between Euler or Runge-Kutta integration method. 
-dt is used to choose integration increment. 
While this model is able to generate spiking activity at high step values (0.01), when temporal study is required, it is recommended to use Euler at least at 0.001 or Runge-Kutta method for more precission results. 

### Example run complete circuit (no ramp)
For an example simulation of 10 seconds you can run the model with the following arguments:

	./feeding_cpg -connection 3 -file_name ./data/complete -integrator -e -dt 0.001 -c_so 8.5 -c_n1m 6 -c_n2v 2 -c_n3t 0 -secs_dur 10

This will run a simulation with the complete model (connection=3) integrated by Euler method using 0.001 as time step for 10 seconds. Each neuron will be stimulated by 8.5, 6, 2 and 0 values.

To see the result run 

	make plot_last or sh ./utils/plot_last.sh
	
	
(see section Plot Utils for more options)

### Example run with ramp
For an example of neuron stimulation by a ramp current injection the current value for the neuron which	you want to stimulate must be -1. 

	./feeding_cpg -connection 3 -file_name ./data/n1m-driven -integrator -e -dt 0.001 -c_so 8.5 -c_n1m -1 -c_n2v 2 -c_n3t 0 -stim_dur 4.6 -stim_inc 0.5 -MIN_c 0 -MAX_c 10 -rounds 4

This will run a simulation with the complete model (connection=3) integrated by Euler method using 0.001 as time step, where neuron N1M is stimulated by the ramp value, which will go from 0 (MIN_c) to 10 (MAX_c) and back twice (rounds=4), the duration of the stimulous is 4.6 secs (stim_dur) and the increment of the current value is 0.5 (stim_inc). Since -secs_dur is not given, the duration time depends on the ramp, when given, the time will be fixed to the specified value.

Note: More than one neuron can be stimulated by the ramp but it is usually not convinient for the rhythm. The ideal is that only one neuron has the stimulation by ramp in each simulation, whereas the rest of them are stimulated by a constant current value. 

To see the result run 

	make plot_last or sh ./utils/plot_last.sh
	
(see section Plot Utils for more options)

### Example run simulating satiated behaviour
For an example simulation of 20 seconds with satiated situation in the feeding CPG simulation between 7s (satiated_ini) and 13s (satiated_end) you can run the model with the following arguments:

	./feeding_cpg -connection 3 -file_name ./data/feeding -integrator -e -dt 0.001 -c_so 8.5 -c_n1m 6 -c_n2v 2 -c_n3t 0 -secs_dur 20 -satiated_ini 7 -satiated_end 13.

After finishing the satiated feeding mode and starting over again the rhythm, current values will stay as they were before the begining of this scenario.

To see the result run 

	make plot_last or sh ./utils/plot_last.sh

(see section Plot Utils for more options)

#### Connectivity options
##### N1M-N2v

	./feeding_cpg -connection 1 -file_name ./data/complete -integrator -e -dt 0.001 -c_n1m 6 -c_n2v 2 -secs_dur 10

Going back to the first example, we can use the same configuration in order to run a simulation where only N1M and N2v are connected (connection=1). The file generated only contains the V value of these two neurons. 

##### N1M-N2v-N3t

	./feeding_cpg -connection 3 -file_name ./data/complete -integrator -e -dt 0.001 -c_n1m 6 -c_n2v 2 -c_n3t 0 -secs_dur 10

Changing "-connection" argument and adding -c_n3t value, N1M, N2v and N3t will be connected during the simulation. The file generated only contains the V value of these three neurons. 


##### Complete circuit: N1M-N2v-N3t-SO 

	./feeding_cpg -connection 3 -file_name ./data/complete -integrator -e -dt 0.001 -c_so 8.5 -c_n1m 6 -c_n2v 2 -c_n3t 0 -secs_dur 10

See Example run complete circuit (no ramp).

##### Complete circuit with synapses
	./feeding_cpg -connection 4 -file_name ./data/complete_syn -integrator -e -dt 0.001 -c_so 8.5 -c_n1m 6 -c_n2v 2 -c_n3t 0 -secs_dur 10

The recorded file includes each neuron voltage in the soma and the synaptic current received (connection=4). 

##### Complete circuit without ramp
	./feeding_cpg -connection 5 -file_name ./data/complete_syn -integrator -e -dt 0.001 -c_so 8.5 -c_n1m 6 -c_n2v 2 -c_n3t 0 -secs_dur 10

The recorded file includes each neuron voltage and ramp stimulation current is not included (connection>=5). 

	
### Plot Utils 
In directory utils you can find some code in python to visualized the generated data during the simulation. 
	
Plot last file generated in data/ using Python3.6. It is plotted in both modes: no-spike and spike. 

	make plot_last or sh ./utils/plot_last.sh

Plot file with file_name. If path is not specified, file generated in data/ using Python3.6 will be plotted.

	python3 plot.py <path> <file_name> 

Plot file with file_name in path with the spikes detected. If path is not specified, file generated in data/ using Python3.6 will be plotted.

	python3 plot_with_spikes.py <path> <file_name>

Warning:
This plotting tool might have difficulties plotting large files, you can use other plotting tools such as gnunplot.

### Developed by ###
Alicia Garrido Peña, PhD in Universidad Autónoma de Madrid. Grupo de Neurocomputación biológica (GNB)

http://www.ii.uam.es/~gnb/

Special thanks to the rest of the members in the group for their help and advice.

### Contact ###

For any questions, please, contact me on:

alicia.garrido@uam.es


