# Developed by Alicia Garrido Peña (2020)
#
# Plotting tools for Lymnaea CPG Simulator Model. 
#
# Implementation of the Lymnaea feeding CPG originally proposed by Vavoulis et al. (2007). Dynamic control of a central pattern generator circuit: A computational model of the snail feeding network. European Journal of Neuroscience, 25(9), 2805–2818. https://doi.org/10.1111/j.1460-9568.2007.05517.x
# and used in study of dynamical invaraiants in Alicia Garrido-Peña, Irene Elices and Pablo Varona (2020). Characterization of interval variability in the sequential activity of a central pattern generator model. Neurocomputing 2020.
#
# Please, if you use this implementation cite the two papers above in your work. 
############################################################################################


import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 50})

import sys
import os
import pandas as pd

if len(sys.argv) >2:
	path = sys.argv[1]
	file_name = sys.argv[2]
else:
	print("Error: No file specified \n Format: <path> <file_name>")
	exit()

print("Ploting file from ",file_name)

path = path+file_name


f = open(path)
headers = f.readline().split()
print(headers)
f.close()

data = pd.read_csv(path, delimiter = " ", names=headers,skiprows=1,low_memory=False)

rows = data.shape[1]

colors = ['teal', 'brown', 'blue', 'green','maroon','teal', 'brown', 'blue', 'green','maroon']

plt.figure(figsize=(30,30))
for i in range(1,rows):
	if(i==1):
		ax1 = plt.subplot(rows,1,i)
	else:	
		plt.subplot(rows,1,i,sharex=ax1)

	if(i==rows-1):
		plt.xlabel("Time (ms)")

	if(headers[i]=='c'):
		plt.ylabel("Current", multialignment='center')	
	else:
		plt.ylabel("Voltage\n(mV)", multialignment='center')
	plt.plot(data['t'],data[headers[i]],color=colors[i-1])

	plt.title(headers[i])

plt.tight_layout()

plt.savefig("./images/"+file_name[:-3]+"eps",format='eps')
plt.show()


