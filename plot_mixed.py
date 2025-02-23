import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot(data1,data2,data3):
    size_img = 1
    plt.rcParams.update({'font.size': 12})
    # plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    alpha = data1[:,0]

    n1 = data1[:,2]
    n2 = data2[:,2]
    n3 = data3[:,2]

    t1 = data1[:,3]
    t2 = data2[:,3]
    t3 = data3[:,3]

    ball_tickness = 2

    # Plot alpha vs n
    plt.plot(alpha,n1,":o",ms=ball_tickness,lw=0.5,label="$\epsilon: 0.05$")
    plt.plot(alpha,n2,":o",ms=ball_tickness,lw=0.5,label="$\epsilon: 0.005$")
    plt.plot(alpha,n3,":o",ms=ball_tickness,lw=0.5,label="$\epsilon: 0.0005$")
    plt.xlabel("$\\alpha$",fontsize=16)
    plt.ylabel("$n$",fontsize=16)
    plt.xticks(np.arange(0.5,1.01,0.1))
    plt.yticks(np.arange(0,601,100))
    plt.ylim(-50,600)
    plt.tick_params(axis='both',direction='in')
    plt.legend(fontsize=11,loc="lower left")
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.3)
    plt.savefig("alpha_vs_n.pdf",bbox_inches = "tight")
    plt.close()

    # Plot alpha vs t
    plt.plot(alpha,t1,":o",ms=ball_tickness,lw=0.5,label="$\epsilon: 0.05$")
    plt.plot(alpha,t2,":o",ms=ball_tickness,lw=0.5,label="$\epsilon: 0.005$")
    plt.plot(alpha,t3,":o",ms=ball_tickness,lw=0.5,label="$\epsilon: 0.0005$")
    plt.xlabel("$\\alpha$",fontsize=16)
    plt.ylabel("$t$",fontsize=16)
    plt.xticks(np.arange(0.5,1.01,0.1))
    plt.yticks(np.arange(0,3601,900))
    plt.ylim(-200,4000)
    plt.tick_params(axis='both',direction='in')
    plt.legend(fontsize=11,loc="lower left")
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.3)
    plt.savefig("alpha_vs_t.pdf",bbox_inches = "tight")


    # plt.xticks(np.arange(0.5,1.01,0.1))
    # plt.yticks(np.arange(0,601,100))
    # plt.xlabel("$\\alpha$",fontsize=16)
    # plt.ylabel("$n$",fontsize=16)
    # plt.tick_params(axis='both',direction='in')
    # plt.savefig("alpha_vs_n.pdf",bbox_inches = "tight")
    # plt.close()

    # # Plot alpha vs t
    # plt.plot(alpha,t,":o",ms=3,lw=1)
    # plt.xticks(np.arange(0.5,1.01,0.1))
    # plt.yticks(np.arange(0,3601,900))
    # plt.ylim(-200,4000)
    # plt.xlabel("$\\alpha$",fontsize=16)
    # plt.ylabel("$t$",fontsize=16)
    # plt.savefig(t_plot[op-1],bbox_inches = "tight")
    # plt.close()


folder = "test_alpha_0.5_1"

data1 = np.loadtxt(folder+"/output_tol_0.05.txt",delimiter=" ")
data2 = np.loadtxt(folder+"/output_tol_0.005.txt",delimiter=" ")
data3 = np.loadtxt(folder+"/output_tol_0.0005.txt",delimiter=" ")

plot(data1,data2,data3)

