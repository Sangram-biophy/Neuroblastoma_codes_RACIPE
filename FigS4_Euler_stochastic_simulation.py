import numpy as np
import matplotlib.pyplot as plt
import os
import math as m
from scipy.integrate import odeint
import pandas as pd
from matplotlib.pyplot import figure


def hill(X,X0,lamb,n):
    return (1/lamb)*((X0**n)/(X0**n+X**n)) + ((X**n)/(X0**n+X**n))
def hilld(X,X0,lamb,n):
    return (X0**n)/(X0**n+X**n) + lamb*((X**n)/(X0**n+X**n))


def odes(x):
    H_PHOX2B_SOX10=hilld(x[0],x_PHOX2B_SOX10,l_PHOX2B_SOX10,n_PHOX2B_SOX10);
    H_SOX10_PHOX2B=hill(x[1],x_SOX10_PHOX2B,l_SOX10_PHOX2B,n_SOX10_PHOX2B);
    H_PRRX1_PHOX2B=hilld(x[2],x_PRRX1_PHOX2B,l_PRRX1_PHOX2B,n_PRRX1_PHOX2B);
    H_PRRX1_SOX10=hill(x[2],x_PRRX1_SOX10,l_PRRX1_SOX10,n_PRRX1_SOX10);
    H_SOX10_PRRX1=hill(x[1],x_SOX10_PRRX1,l_SOX10_PRRX1,n_SOX10_PRRX1);
    H_SOX10_SOX10=hill(x[1],x_SOX10_SOX10,l_SOX10_SOX10,n_SOX10_SOX10);
    H_PHOX2B_PHOX2B=hill(x[0],x_PHOX2B_PHOX2B,l_PHOX2B_PHOX2B,n_PHOX2B_PHOX2B);
    H_GATA3_PHOX2B=hill(x[3],x_GATA3_PHOX2B,l_GATA3_PHOX2B,n_GATA3_PHOX2B);
    H_PHOX2B_GATA3=hill(x[0],x_PHOX2B_GATA3,l_PHOX2B_GATA3,n_PHOX2B_GATA3);
    H_GATA3_GATA3=hill(x[3],x_GATA3_GATA3,l_GATA3_GATA3,n_GATA3_GATA3);
    H_PRRX1_GATA3=hilld(x[2],x_PRRX1_GATA3,l_PRRX1_GATA3,n_PRRX1_GATA3);
    H_SNAI2_SOX10=hill(x[4],x_SNAI2_SOX10,l_SNAI2_SOX10,n_SNAI2_SOX10);
    H_SOX10_SNAI2=hill(x[1],x_SOX10_SNAI2,l_SOX10_SNAI2,n_SOX10_SNAI2);
    H_SNAI2_PRRX1=hill(x[4],x_SNAI2_PRRX1,l_SNAI2_PRRX1,n_SNAI2_PRRX1);
    H_PRRX1_SNAI2=hill(x[2],x_PRRX1_SNAI2,l_PRRX1_SNAI2,n_PRRX1_SNAI2);
    H_SNAI2_SNAI2=hill(x[4],x_SNAI2_SNAI2,l_SNAI2_SNAI2,n_SNAI2_SNAI2);
    H_HAND2_PHOX2B=hill(x[5],x_HAND2_PHOX2B,l_HAND2_PHOX2B,n_HAND2_PHOX2B);
    H_PHOX2B_HAND2=hill(x[0],x_PHOX2B_HAND2,l_PHOX2B_HAND2,n_PHOX2B_HAND2);
    H_HAND2_GATA3=hill(x[5],x_HAND2_GATA3,l_HAND2_GATA3,n_HAND2_GATA3);
    H_GATA3_HAND2=hill(x[3],x_GATA3_HAND2,l_GATA3_HAND2,n_GATA3_HAND2);
    H_HAND2_HAND2=hill(x[5],x_HAND2_HAND2,l_HAND2_HAND2,n_HAND2_HAND2);


    dPHOX2B = g_PHOX2B*H_SOX10_PHOX2B*H_PRRX1_PHOX2B*H_PHOX2B_PHOX2B*H_GATA3_PHOX2B*H_HAND2_PHOX2B - k_PHOX2B*x[0];
    dSOX10 = g_SOX10*H_PHOX2B_SOX10*H_PRRX1_SOX10*H_SOX10_SOX10*H_SNAI2_SOX10 - k_SOX10*x[1];
    dPRRX1 = g_PRRX1*H_SOX10_PRRX1*H_SNAI2_PRRX1 - k_PRRX1*x[2];
    dGATA3 = g_GATA3*H_PHOX2B_GATA3*H_GATA3_GATA3*H_PRRX1_GATA3*H_HAND2_GATA3 - k_GATA3*x[3];
    dSNAI2 = g_SNAI2*H_SOX10_SNAI2*H_PRRX1_SNAI2*H_SNAI2_SNAI2 - k_SNAI2*x[4];
    dHAND2 = g_HAND2*H_PHOX2B_HAND2*H_GATA3_HAND2*H_HAND2_HAND2 - k_HAND2*x[5];
    
    return [dPHOX2B, dSOX10, dPRRX1, dGATA3, dSNAI2, dHAND2]


g_PHOX2B = 378.334864;
g_SOX10 = 55.753313;
g_PRRX1 = 70.443192;
g_GATA3 = 98.642177;
g_SNAI2 = 78.89362;
g_HAND2 = 4.639991;


k_PHOX2B = 0.161111;
k_SOX10 = 0.498492;
k_PRRX1 = 0.38904;
#k_PRRX1=2.25;
k_GATA3 = 0.28408;
k_SNAI2 = 0.506437;
k_HAND2 = 0.522097;


x_PHOX2B_SOX10 = 6.437256;
x_SOX10_PHOX2B = 1.166929;
x_PRRX1_PHOX2B = 20.919588;
x_PRRX1_SOX10 = 1.980271;
x_SOX10_PRRX1 = 1.42407;
x_SOX10_SOX10 = 1.587647;
x_PHOX2B_PHOX2B = 3.373513;
x_GATA3_PHOX2B = 1.9518;
x_PHOX2B_GATA3 = 6.229318;
x_GATA3_GATA3 = 0.124213;
x_PRRX1_GATA3 = 8.07046;
x_SNAI2_SOX10 = 0.084409;
x_SOX10_SNAI2 = 2.039205;
x_SNAI2_PRRX1 = 4.546443;
x_PRRX1_SNAI2 = 4.725253;
x_SNAI2_SNAI2 = 0.971109;
x_HAND2_PHOX2B = 2.59919;
x_PHOX2B_HAND2 = 1.075751;
x_HAND2_GATA3 = 4.635078;
x_GATA3_HAND2 = 1.290636;
x_HAND2_HAND2 = 1.937974;


l_PHOX2B_SOX10 = 0.010113;
l_SOX10_PHOX2B = 49.674949;
l_PRRX1_PHOX2B = 0.013027;
l_PRRX1_SOX10 = 70.651465;
l_SOX10_PRRX1 = 65.988884;
l_SOX10_SOX10 = 10.644229;
l_PHOX2B_PHOX2B = 98.135254;
l_GATA3_PHOX2B = 89.640035;
l_PHOX2B_GATA3 = 65.390146;
l_GATA3_GATA3 = 39.656928;
l_PRRX1_GATA3 = 0.023987;
l_SNAI2_SOX10 = 40.051565;
l_SOX10_SNAI2 = 23.532513;
l_SNAI2_PRRX1 = 55.037178;
l_PRRX1_SNAI2 = 98.171792;
l_SNAI2_SNAI2 = 87.974437;
l_HAND2_PHOX2B = 77.957324;
l_PHOX2B_HAND2 = 27.563915;
l_HAND2_GATA3 = 97.963846;
l_GATA3_HAND2 = 68.971081;
l_HAND2_HAND2 = 38.348446;


n_PHOX2B_SOX10 = 1;
n_SOX10_PHOX2B = 5;
n_PRRX1_PHOX2B = 6;
n_PRRX1_SOX10 = 6;
n_SOX10_PRRX1 = 6;
n_SOX10_SOX10 = 1;
n_PHOX2B_PHOX2B = 5;
n_GATA3_PHOX2B = 1;
n_PHOX2B_GATA3 = 5;
n_GATA3_GATA3 = 4;
n_PRRX1_GATA3 = 6;
n_SNAI2_SOX10 = 4;
n_SOX10_SNAI2 = 5;
n_SNAI2_PRRX1 = 4;
n_PRRX1_SNAI2 = 3;
n_SNAI2_SNAI2 = 3;
n_HAND2_PHOX2B = 5;
n_PHOX2B_HAND2 = 2;
n_HAND2_GATA3 = 1;
n_GATA3_HAND2 = 2;
n_HAND2_HAND2 = 4;

M=10;


time_dyn=np.reshape([0,0,0,0,0,0,0],(1,7))
dt=0.01
t_tot=200
noise=5
nn=int(t_tot/dt)+1
#time_dyn=np.zeros([nn,3],dtype=float)
t=0
i=0
X0=np.array([(g_PHOX2B/k_PHOX2B),(g_SOX10/k_SOX10),(g_PRRX1/k_PRRX1),g_GATA3/k_GATA3,g_SNAI2/k_SNAI2,g_HAND2/k_HAND2])
#X0=np.array([10,10,10,10,10,10])
while t<=t_tot:
    X=X0+np.array(odes(X0))*dt+m.sqrt(noise)*np.random.normal(loc=0.0, scale=1.0, size=(6))
    #print(X)
   # time_dyn[i,:]=[t,X[0],X[1]]
    if np.any(X<0.0):
        neg=np.where(X<0.0)[0]
        for ne in range(len(neg)):
            X[neg[ne]]=0.0
        #time_dyn=np.append(time_dyn,np.reshape([t,0,0,0,0,0,0],(1,7)),axis=0)
    
    time_dyn=np.append(time_dyn,np.reshape([t,X[0],X[1],X[2],X[3],X[4],X[5]],(1,7)),axis=0)

    X0=X
    t=t+dt
    #i=i+1
time_dyn=np.delete(time_dyn,0,0)


df = pd.DataFrame(time_dyn, columns = ['t','PHOX2B', 'SOX10', 'PRRX1', 'GATA3', 'SNAI2', 'HAND2'])


figure(figsize=(15, 7.5), dpi=100)
plt.rcParams.update({'font.size': 30})
plt.plot(df["t"],df["PHOX2B"],label="PHOX2B")
plt.plot(df["t"],df["SOX10"],label="SOX10")
plt.plot(df["t"],df["PRRX1"],label="PRRX1")
plt.plot(df["t"],df["GATA3"],label="GATA3")
plt.plot(df["t"],df["SNAI2"],label="SNAI2")
plt.plot(df["t"],df["HAND2"],label="HAND2")
plt.rcParams.update({'font.size': 30})
plt.xlim(0,200)
plt.ylim(0,300)
plt.legend(bbox_to_anchor =(1, 1))
plt.show()


figure(figsize=(15, 7.5), dpi=100)
plt.rcParams.update({'font.size': 30})
xx="SOX10"
yy="PHOX2B"
plt.plot(df["t"],df[xx],label=xx)
plt.plot(df["t"],df[yy],label=yy)
plt.rcParams.update({'font.size': 30})
plt.xlim(175,200)
plt.ylim(0,120)
plt.legend(bbox_to_anchor =(1, 1))
plt.show()

df.to_csv('euler_stochastic.csv', index=False)







