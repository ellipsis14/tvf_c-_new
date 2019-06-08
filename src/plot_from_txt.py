import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data1 = pd.read_csv('convergence.txt', sep=',', header = None)
data2 = pd.read_csv('logFile.txt', sep=',', header = None)


data1.columns = ['l12', 'l13', 'h12', 'h13', 'dg12', 'dg13']
data2.columns = ['t', 'it', 'nt', 'o', 'b', '']


# Find the ratio

data1['l2'] = data1.l12 /data1.l13 
data1['h1'] = data1.h12 /data1.h13 
data1['COIPDG'] = data1.dg12 /data1.dg13 
data2['t'] = data2.t

#

fig, ax1 = plt.subplots()
ax1.plot(data2.t, data1.l2, color='red')
ax2 = ax1 
ax2.plot(data2.t, data1.h1, color='blue')
ax3 = ax1
ax3.plot(data2.t, np.log2(data1.COIPDG), color='green')
plt.title('Convergence_plot')
plt.legend()
plt.show()
