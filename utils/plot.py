import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 50})

import sys
import os
import pandas as pd

if len(sys.argv) >2:
	path = sys.argv[2]
else:
	path = "./data/"
if len(sys.argv) >1:
	file_name = sys.argv[1]
else:
	print("Error: No file specified")
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


