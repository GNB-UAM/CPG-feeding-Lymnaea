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
	file_name = "prueba.asc"

print("Ploting file from ",file_name)

path = path+file_name
# path = file_name


f = open(path)
headers = f.readline().split()
print(headers)
f.close()

# data = np.loadtxt(path,skiprows=1)



data = pd.read_csv(path, delimiter = " ", names=headers,skiprows=1,low_memory=False)


# print("Read!")
rows = data.shape[1]
# rows = len(headers)

colors = ['teal', 'brown', 'blue', 'green','maroon','teal', 'brown', 'blue', 'green','maroon']
# colors = ['teal', 'teal','brown','brown', 'blue','blue', 'green', 'green','maroon','maroon']

plt.figure(figsize=(30,30))
for i in range(1,rows):
	# print(i)
	if(i==1):
		ax1 = plt.subplot(rows,1,i)
	else:	
		plt.subplot(rows,1,i,sharex=ax1)

	# if(i%2==1):
	# 	# header= headers[i]+"\n(mV)"
	# 	plt.ylabel("Voltage\n(mV)", multialignment='center')
	# else:
	# 	# header= headers[i]+"\n(nA)"
	# 	plt.ylabel("Current\n(nA)", multialignment='center')
	# plt.plot(data[1:,0],data[1:,i])
	# plt.xlabel("Iterations")
	if(i==rows-1):
		plt.xlabel("Time (ms)")
	# if(i== rows//2):
	# 	plt.ylabel("Voltage (mV)", multialignment='center')	

	if(headers[i]=='c'):
		plt.ylabel("Current", multialignment='center')	
	else:
		plt.ylabel("Voltage\n(mV)", multialignment='center')
	# plt.plot(data['t'],data[headers[i]],color=colors[i%len(colors)])
	plt.plot(data['t'],data[headers[i]],color=colors[i-1])
	# plt.plot(data[headers[i]],color=colors[i-1])
	plt.title(headers[i])
	# plt.title(header)


# plt.ylabel("Voltage (mV)")
# plt.text(-0.2, 0.7, 'Voltage (mV)', ha='center', va='center', rotation='vertical')

# fig, axes = plt.subplots(nrows=rows, ncols=1)
# fig.tight_layout()
# plt.subplots_adjust(hspace = .001)
plt.tight_layout()
# plt.savefig("./images/"+file_name[:-3]+"png")

plt.savefig("./images/"+file_name[:-3]+"eps",format='eps')
plt.show()




# for i in range(1,rows):
# 	print(i)
# 	plt.figure(figsize=(15,5))
# 	# plt.xlabel("Iterations")
# 	plt.xlabel("Time (ms)")
# 	plt.ylabel("Voltage \n(mV)")
# 	plt.plot(data['t'],data[headers[i]])
# 	# plt.xlim(0,8000)
# 	# plt.title(headers[i])
# 	# plt.title("dt = 0.001")
# 	plt.tight_layout()

# 	# plt.savefig("./images/"+file_name[:-4]+headers[i]+"_activity"+".png")
# 	# plt.show()
	
# os.system("rm "+path)

# for i in range(1,4):
# 	plt.plot(data[1:,0],data[1:,i],label=headers[i])
# 	plt.legend()

# plt.show()



# plt.subplot(2,1,1)
# plt.plot(data[1:,0],data[1:,1])
# plt.title("Irregular")
# plt.subplot(2,1,2)
# plt.plot(data[1:,0],data[1:,3])
# plt.title("Regular")
# plt.show()