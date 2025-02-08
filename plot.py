import numpy as np
import matplotlib.pyplot as plt

size_img = 0.6
plt.rcParams.update({'font.size': 12})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(axis='both',direction='in')

data1 = np.loadtxt("alpha05/output_1_10.txt",delimiter=" ")
data2 = np.loadtxt("alpha05/output_10_1000.txt",delimiter=" ")
data3 = np.loadtxt("alpha05/output_1000_100000.txt",delimiter=" ")

x1 = data1[:,0]
y1 = data1[:,2]

x2 = data2[:,0]
y2 = data2[:,2]

x3 = data3[:,0]
y3 = data3[:,2]

plt.plot(x1,y1,":o",ms=3,lw=1)
plt.xlabel("$n$",fontsize=16)
plt.ylabel("$s_n$",fontsize=16)
plt.savefig("alpha05/mnqfig5a.pdf",bbox_inches = "tight")
plt.close()

plt.plot(x2,y2,":o",ms=3,lw=1)
plt.xlabel("$n$",fontsize=16)
plt.ylabel("$s_n$",fontsize=16)
plt.savefig("alpha05/mnqfig5b.pdf",bbox_inches = "tight")
plt.close()

plt.plot(x3,y3,":o",ms=3,lw=1)
plt.xlabel("$n$",fontsize=16)
plt.ylabel("$s_n$",fontsize=16)
plt.savefig("alpha05/mnqfig5c.pdf",bbox_inches = "tight")
