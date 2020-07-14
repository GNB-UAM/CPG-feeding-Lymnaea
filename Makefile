LIBDIR=./include
SRCDIR=./src/
PLOTDIR=./utils/

LDIR=-L$(LIBDIR)
IDIR=-I$(LIBDIR)
CFLAGS=-Wall
COPT=-O2
CC=g++ -std=c++11

all: simulation


simulation: $(SRCDIR)lymnaea_main.cpp $(SRCDIR)vavoulis_neuron.cpp $(SRCDIR)vavoulis_synapse.cpp $(SRCDIR)ramp_generator.cpp $(SRCDIR)cpg_simulator.cpp
	$(CC) $(CFLAGS) $(COPT) $(SRCDIR)lymnaea_main.cpp $(SRCDIR)vavoulis_neuron.cpp $(SRCDIR)vavoulis_synapse.cpp $(SRCDIR)ramp_generator.cpp $(SRCDIR)cpg_simulator.cpp -o feeding_cpg -lm -I$(LIBDIR)

run_default: simulation 
	./feeding_cpg -connection 3 -file_name ./data/complete -integrator -e -dt 0.001 -c_so 8.5 -c_n1m 6 -c_n2v 2 -c_n3t 0 -secs_dur 10

plot_last:
	sh $(PLOTDIR)plot_last.sh

config:
	@mkdir -p images data

doxygen:
	doxygen Doxyfile

clean:
	rm -f feeding_cpg *.o 
	rm -f -r html/* latex/*
	rmdir html latex
