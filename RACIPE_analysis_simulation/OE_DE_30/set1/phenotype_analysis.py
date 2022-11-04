#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats 


# In[2]:


def stable_states(path,file_start):
    num_stable_state=[]
    file=[]
    for filename in os.listdir(path):
        if filename[-14:-11]=="sol":
            num_stable_state.append(filename.split("_")[4].split(".")[0])
            expr=filename.split("_")[1]+filename.split("_")[2]
            file.append(filename)
#         print(q)
    return(num_stable_state,file,expr)


# In[3]:


def nodes(path):
    node_name=[]
    expression={}
    for filename in os.listdir(path):
        if filename.endswith(".cfg"):
            with open(path+filename,'r') as g:
                flag=-1
                for line in g:
                    a=line.split()
                    if line.startswith("number_of_"):
                        a=line.split()
                        expression[a[0]]=int(a[1])
                    if a[0]=="NumberOfGenes":
                        node=int(a[1])
                        flag=0
                        continue
                    if flag >=0 and flag<node:
                        flag+=1
                        node_name.append(a[1])
    for exprr, number in expression.items():
        if number==1:
            expression_type=exprr.split("_")[2][0:2]
#expression_type
    return node,node_name,expression_type


# In[4]:


def fold_change(path,expression_type):
    for filename in os.listdir(path):
            if filename.endswith(".cfg"):
                with open(path+filename,'r') as g:
                    for line in g:
                        if line.startswith(expression_type+"_ID"):
                            expr_node=node_name[int(line.split()[1])-1]
                        if line.startswith(expression_type+"_Fold_Change"):
                            expr_fch=int(float(line.split()[1]))
    return expr_node, expr_fch


# In[5]:


def state_count(state, file):
    for i in range(len(state)):
        with open(path+file[i],'r') as f:
                count1=0
                for line in f:
                    count1=count1+1
        print("state=",state[i],count1)


# In[ ]:





# In[6]:


def collated_data(nodes):
    data=np.empty((1,nodes),dtype=float)
    for i in range(len(state)):
        with open(path+file[i],'r') as f:
            for line in f:
                for j in range(int(state[i])):
#             print(line.split("\t")[nodes*i+2:nodes*(i+1)+2])
                    q=line.strip().split("\t")[nodes*j+2:nodes*(j+1)+2]
                    q=np.array([float(k) for k in q])
            #q=np.array(q)
                    data=np.append(data,np.reshape(q,(1,nodes)),axis=0)
                #print(q)
    data=np.delete(data,0,0)
    return data


# In[7]:


def hist_kde_plot(histpath,z_score):
    #histpath=path+"\\analysis\\hist_kde\\"
    if not os.path.exists(histpath):
        os.mkdir(histpath)
    #os.mkdir(histpath)
    fig, ax = plt.subplots(figsize = (15, 7.5),constrained_layout=True)
    ax.minorticks_on()
    for i in range(z_score.shape[1]):
        sns.distplot(z_score[z_score.columns[i]],color='purple')
        plt.ylim(0, 1)
        plt.rcParams.update({'font.size': 30})
    #plt.legend(prop={"size":20})
        plt.savefig(histpath+z_score.columns[i]+"_kde_hist"+'.png',bbox_inches="tight",dpi=600)
        plt.clf()


# In[8]:


def scatter_color(x,y,z,point_size,cbar_label,xlab,ylab):
    fig, ax = plt.subplots(figsize = (15, 7.5),constrained_layout=True)
    plt.scatter(x,y,s=point_size, edgecolor='black', linewidth=0, c=z,cmap=plt.cm.get_cmap("bwr",20))
    cbar=plt.colorbar(orientation="vertical")
    cbar.set_label(label=cbar_label,size=25)
    plt.clim(-3,3)
    cbar.ax.tick_params(labelsize=20)
    plt.rcParams.update({'font.size': 30})
    plt.xlabel(xlab)
    plt.ylabel(ylab)


# In[9]:


path=os.getcwd()+"/"
file_start="core_OE"
state,file, expr=stable_states(path,file_start)
node, node_name,expression_type=nodes(path)
expr_node, expr_fch=fold_change(path,expression_type)


# In[10]:


collated=collated_data(node)
head=""
for i in range(len(node_name)):
    head=head + "\t" + node_name[i]
head=head.strip()
collated_file_name='collated_core'+expr
np.savetxt(collated_file_name+'.txt',collated,fmt='%1.6f', delimiter='\t', header=head, comments='')


# In[11]:


z_score=pd.DataFrame(stats.zscore(collated,axis=0),columns=node_name)
core_path='/home/user/Documents/Neuroblastoma_MKJlab/trial/CORE/'
core_data=pd.read_csv(core_path+"collated_core.txt",sep="\t")
collated_OE=pd.DataFrame(collated,columns=node_name)
OE_diff=(collated_OE-core_data.mean())/core_data.std()


# In[12]:


analysis_folder=os.getcwd()+"/analysis/"
if not os.path.exists(analysis_folder):
    os.mkdir(analysis_folder)


# In[13]:


histpath=path+"kde_plot/"
hist_kde_plot(histpath,z_score)


# In[14]:


path=os.getcwd()+"/analysis/"
y_axis="PRRX1"
x_axis="PHOX2B"
z_axis="SOX10"
scatter_color(OE_diff[x_axis],OE_diff[y_axis],OE_diff[z_axis],point_size=7,cbar_label=z_axis,xlab=x_axis,ylab=y_axis)
plt.clim(-2,2)
plt.title(expression_type+"_"+expr_node+"_"+str(expr_fch)+"fold")
plt.savefig(path+x_axis+"_"+y_axis+"_"+z_axis+"_"+expression_type+"_"+expr_node+"_"+str(expr_fch)+"_fold.png",bbox_inches="tight",dpi=600)


# In[15]:


ll=0
hl=0
lh=0
hh=0
for i in range(len(OE_diff)):
    if OE_diff[x_axis][i]>0.1 and OE_diff[y_axis][i]<0:
        hl=hl+1
    elif OE_diff[x_axis][i]<0.1 and OE_diff[y_axis][i]<0:
        ll=ll+1
    elif OE_diff[x_axis][i]<0.1 and OE_diff[y_axis][i]>0:
        lh=lh+1
    elif OE_diff[x_axis][i]>0.1 and OE_diff[y_axis][i]>0:
        hh=hh+1
phenotype=['hl','ll','lh','hh']
per_pheno=np.array([hl,ll,lh,hh])/len(OE_diff)


# In[16]:


fig, ax = plt.subplots(figsize = (15, 7.5),constrained_layout=True)
plt.bar(phenotype,per_pheno)
ax.set_ylabel('Frequency')
ax.set_xlabel('Phenotype')
plt.title(expression_type+"_"+expr_node+"_"+str(expr_fch)+"fold")
plt.savefig(path+expression_type+"_"+expr_node+"_"+str(expr_fch)+"_hist_phenotype.png",bbox_inches="tight",dpi=600)


# In[ ]:




