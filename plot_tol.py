import numpy as np
import matplotlib.pyplot as plt

def plot(data,op):
    alpha = data[:,0]
    n = data[:,2]
    t = data[:,3]

    n_plot = ["alpha_vs_n_0.005.pdf","alpha_vs_n_0.0005.pdf"]

    t_plot = ["alpha_vs_t_0.005.pdf","alpha_vs_t_0.0005.pdf"]


    # Plot alpha vs n
    plt.plot(alpha,n,":o",ms=3,lw=1)
    plt.xticks(np.arange(0.5,1.01,0.1))
    plt.yticks(np.arange(0,601,100))
    plt.xlabel("$\\alpha$",fontsize=16)
    plt.ylabel("$n$",fontsize=16)
    plt.tick_params(axis='both',direction='in')
    plt.savefig(n_plot[op-1],bbox_inches = "tight")
    plt.close()

    # Plot alpha vs t
    plt.plot(alpha,t,":o",ms=3,lw=1)
    plt.xticks(np.arange(0.5,1.01,0.1))
    plt.yticks(np.arange(0,3601,900))
    plt.ylim(-200,4000)
    plt.xlabel("$\\alpha$",fontsize=16)
    plt.ylabel("$t$",fontsize=16)
    plt.tick_params(axis='both',direction='in')
    plt.savefig(t_plot[op-1],bbox_inches = "tight")
    plt.close()

    

size_img = 0.6
plt.rcParams.update({'font.size': 12})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(axis='both',direction='in')

folder = "test_alpha_0.5_1"

data1 = np.loadtxt(folder+"/output_tol_0.005.txt",delimiter=" ")
data2 = np.loadtxt(folder+"/output_tol_0.0005.txt",delimiter=" ")

plot(data1,1)
plot(data2,2)
