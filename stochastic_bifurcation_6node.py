import numpy as np
import matplotlib.pyplot as plt
import os
import math as m
from scipy.integrate import odeint
import pandas as pd
import seaborn as sns

def hill(X,X0,lamb,n):
    return (1/lamb)*((X0**n)/(X0**n+X**n)) + ((X**n)/(X0**n+X**n))
def hilld(X,X0,lamb,n):
    return (X0**n)/(X0**n+X**n) + lamb*((X**n)/(X0**n+X**n))

g_PHOX2B = 378.334864;
g_SOX10 = 55.753313;
g_PRRX1 = 70.443192;
g_GATA3 = 98.642177;
g_SNAI2 = 78.89362;
g_HAND2 = 4.639991;


k_PHOX2B = 0.161111;
k_SOX10 = 0.498492;
#k_PRRX1 = 0.38904;
k_PRRX1=0.1;
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




#init_cond=[(g_PHOX2B/k_PHOX2B),(g_SOX10/k_SOX10),(g_PRRX1/k_PRRX1),g_GATA3/k_GATA3,g_SNAI2/k_SNAI2,g_HAND2/k_HAND2]
init_cond=np.array([10,10,10,10,10,10])
x=init_cond*M

t=0.0
t0=0.0
ti=0.0
tf=20.0
n_point=1000
dt=(tf-ti)/n_point
num=[]
dk=(3.5-0.1)/n_point
#time_dyn=np.zeros((1000,4),dtype=float)
time_dyn=np.empty((1,8),dtype=float)
for i in range(n_point):
    
    kPRRX1=k_PRRX1+i*dk
    #init_cond=[(g_PHOX2B/k_PHOX2B),(g_SOX10/k_SOX10),(g_PRRX1/kPRRX1),g_GATA3/k_GATA3,g_SNAI2/k_SNAI2,g_HAND2/k_HAND2]
#init_cond=np.array([0,98,188,0,158,0])
    x=init_cond*M
    print("starting loop",i)
    t=0
    while t<100:
    
        H_PHOX2B_SOX10=hilld(x[0]/M,x_PHOX2B_SOX10,l_PHOX2B_SOX10,n_PHOX2B_SOX10);
        H_SOX10_PHOX2B=hill(x[1]/M,x_SOX10_PHOX2B,l_SOX10_PHOX2B,n_SOX10_PHOX2B);
        H_PRRX1_PHOX2B=hilld(x[2]/M,x_PRRX1_PHOX2B,l_PRRX1_PHOX2B,n_PRRX1_PHOX2B);
        H_PRRX1_SOX10=hill(x[2]/M,x_PRRX1_SOX10,l_PRRX1_SOX10,n_PRRX1_SOX10);
        H_SOX10_PRRX1=hill(x[1]/M,x_SOX10_PRRX1,l_SOX10_PRRX1,n_SOX10_PRRX1);
        H_SOX10_SOX10=hill(x[1]/M,x_SOX10_SOX10,l_SOX10_SOX10,n_SOX10_SOX10);
        H_PHOX2B_PHOX2B=hill(x[0]/M,x_PHOX2B_PHOX2B,l_PHOX2B_PHOX2B,n_PHOX2B_PHOX2B);
        H_GATA3_PHOX2B=hill(x[3]/M,x_GATA3_PHOX2B,l_GATA3_PHOX2B,n_GATA3_PHOX2B);
        H_PHOX2B_GATA3=hill(x[0]/M,x_PHOX2B_GATA3,l_PHOX2B_GATA3,n_PHOX2B_GATA3);
        H_GATA3_GATA3=hill(x[3]/M,x_GATA3_GATA3,l_GATA3_GATA3,n_GATA3_GATA3);
        H_PRRX1_GATA3=hilld(x[2]/M,x_PRRX1_GATA3,l_PRRX1_GATA3,n_PRRX1_GATA3);
        H_SNAI2_SOX10=hill(x[4]/M,x_SNAI2_SOX10,l_SNAI2_SOX10,n_SNAI2_SOX10);
        H_SOX10_SNAI2=hill(x[1]/M,x_SOX10_SNAI2,l_SOX10_SNAI2,n_SOX10_SNAI2);
        H_SNAI2_PRRX1=hill(x[4]/M,x_SNAI2_PRRX1,l_SNAI2_PRRX1,n_SNAI2_PRRX1);
        H_PRRX1_SNAI2=hill(x[2]/M,x_PRRX1_SNAI2,l_PRRX1_SNAI2,n_PRRX1_SNAI2);
        H_SNAI2_SNAI2=hill(x[4]/M,x_SNAI2_SNAI2,l_SNAI2_SNAI2,n_SNAI2_SNAI2);
        H_HAND2_PHOX2B=hill(x[5]/M,x_HAND2_PHOX2B,l_HAND2_PHOX2B,n_HAND2_PHOX2B);
        H_PHOX2B_HAND2=hill(x[0]/M,x_PHOX2B_HAND2,l_PHOX2B_HAND2,n_PHOX2B_HAND2);
        H_HAND2_GATA3=hill(x[5]/M,x_HAND2_GATA3,l_HAND2_GATA3,n_HAND2_GATA3);
        H_GATA3_HAND2=hill(x[3]/M,x_GATA3_HAND2,l_GATA3_HAND2,n_GATA3_HAND2);
        H_HAND2_HAND2=hill(x[5]/M,x_HAND2_HAND2,l_HAND2_HAND2,n_HAND2_HAND2);
    
    
        propensity=[M*g_PHOX2B*H_SOX10_PHOX2B*H_PRRX1_PHOX2B*H_PHOX2B_PHOX2B*H_GATA3_PHOX2B*H_HAND2_PHOX2B,
                    M*g_SOX10*H_PHOX2B_SOX10*H_PRRX1_SOX10*H_SOX10_SOX10*H_SNAI2_SOX10,
                    M*g_PRRX1*H_SOX10_PRRX1*H_SNAI2_PRRX1,
                    M*g_GATA3*H_PHOX2B_GATA3*H_GATA3_GATA3*H_PRRX1_GATA3*H_HAND2_GATA3,
                    M*g_SNAI2*H_SOX10_SNAI2*H_PRRX1_SNAI2*H_SNAI2_SNAI2,
                    M*g_HAND2*H_PHOX2B_HAND2*H_GATA3_HAND2*H_HAND2_HAND2,
                    k_PHOX2B*x[0],
                    k_SOX10*x[1],
                    kPRRX1*x[2],
                    k_GATA3*x[3],
                    k_SNAI2*x[4],
                    k_HAND2*x[5]]
        a=sum(propensity)
        r=np.random.uniform(size=2)
        tau=(1/a)*(np.log(1/r[0]))
        l=1
        for j in range(len(propensity)):
        #print(sum(propensity[0:j+1]),r[1]*a)
            if sum(propensity[0:j+1])<r[1]*a:
                l=l+1
             #print(l)
            else:
                 break
        num.append(l)
    #print(l)
        if l==1:
            x[0]=x[0]+1
        elif l==2:
            x[1]=x[1]+1
        elif l==3:
            x[2]=x[2]+1
        elif l==4:
            x[3]=x[3]+1
        elif l==5:
            x[4]=x[4]+1
        elif l==6:
            x[5]=x[5]+1
        elif l==7:
            x[0]=x[0]-1
        elif l==8:
            x[1]=x[1]-1
        elif l==9:
            x[2]=x[2]-1
        elif l==10:
            x[3]=x[3]-1
        elif l==11:
            x[4]=x[4]-1
        elif l==12:
            x[5]=x[5]-1
   # print(conc)
        t=t+tau
    #print(t)
        x=np.maximum(0,x)
    time_dyn=np.append(time_dyn,np.reshape([t0,kPRRX1,x[0]/M,x[1]/M,x[2]/M,x[3]/M,x[4]/M,x[5]/M],(1,8)),axis=0)
    t0=t0+dt


df = pd.DataFrame(time_dyn, columns = ['t','k_PRRX1','PHOX2B', 'SOX10', 'PRRX1', 'GATA3', 'SNAI2', 'HAND2'])
df=df.drop([0],axis=0)
df.to_csv('stochastic.csv', index=False)

fig, ax = plt.subplots(figsize = (15, 7.5),constrained_layout=True)
xx="PHOX2B"
yy="GATA3"
ax.plot(df["t"],((df[xx])),color="tab:blue",label=xx)
ax.plot(df["t"],((df[yy])),color="tab:orange",label=yy)
# ax.plot(x,y2)
ax.set_xlabel('Time')
ax.set_ylabel('Concentration')
ax.set_xlim(0,20)

def deg2rad(x):
    return (x/dt)*dk


def rad2deg(x):
    return dt*(x-0)/dk


secax = ax.secondary_xaxis('top', functions=(deg2rad,rad2deg))
secax.set_xlabel('degradation of PRRX1(k_PRRX1)')
secax.set_xlim(0,3.6)
plt.rcParams.update({'font.size': 30})
plt.legend(prop={"size":20})
#plt.savefig('TT_ssa3.png')
plt.ylim(-5,300)
#plt.xlim(0,3.5)
plt.savefig(xx+"_"+yy+'_stochastiic4_plot.png',bbox_inches="tight",dpi=600)
plt.show()

