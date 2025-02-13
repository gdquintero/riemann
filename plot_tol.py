import numpy as np
import matplotlib.pyplot as plt

def plot(data,op):
    alpha = data[:,0]
    n = data[:,2]
    t = data[:,3]

    n_plot = ["alpha_vs_n_0.5.pdf","alpha_vs_n_0.05.pdf","alpha_vs_n_0.005.pdf","alpha_vs_n_0.0005.pdf"]

    t_plot = ["alpha_vs_t_0.5.pdf","alpha_vs_t_0.05.pdf","alpha_vs_t_0.005.pdf","alpha_vs_t_0.0005.pdf"]

    plt.plot(alpha,n,":o",ms=3,lw=1)
    plt.xlabel("$\\alpha$",fontsize=16)
    plt.ylabel("$n$",fontsize=16)
    plt.savefig(n_plot[op-1],bbox_inches = "tight")
    plt.close()
    plt.plot(alpha,t,":o",ms=3,lw=1)
    plt.xlabel("$\\alpha$",fontsize=16)
    plt.ylabel("$t$",fontsize=16)
    plt.savefig(t_plot[op-1],bbox_inches = "tight")
    plt.close()

size_img = 0.6
plt.rcParams.update({'font.size': 12})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(axis='both',direction='in')

data1 = np.loadtxt("output_tol_0.5.txt",delimiter=" ")
data2 = np.loadtxt("output_tol_0.05.txt",delimiter=" ")
data3 = np.loadtxt("output_tol_0.005.txt",delimiter=" ")

plot(data1,1)
plot(data2,2)
plot(data3,3)