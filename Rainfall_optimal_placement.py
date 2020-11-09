#!/usr/bin/env python
# coding: utf-8

# In[3]:


from oct2py import octave


# In[1]:


import shutil
import os
import sys


# In[2]:


os.environ["OCTAVE_EXECUTABLE"] = "C:\\Octave\\Octave-4.4.1\\bin\\octave-cli.exe"


# In[ ]:


octave.addpath('')


# In[4]:


import ax


# In[5]:


from oct2py import Oct2Py
oc = Oct2Py()


# In[6]:


help(oc.fmincon)


# In[11]:


octave.addpath(octave.genpath('C:\\Users\\ramui\\Desktop\\Optimal_Placement\\Code2'))


# In[12]:


octave.fmincon


# In[13]:


help(octave.fmincon)


# In[14]:


help(octave.getmodelmse)


# In[15]:


from scipy.io import loadmat


# In[16]:


yearstart = 1901
yearend = 2000


# In[17]:


monthstart = 6
monthend = 9


# In[18]:


K = 5


# In[19]:


sigmae = 1


# In[20]:


mat_contents = loadmat('raindat0_25_deg.mat')


# In[21]:


hrrainmat= mat_contents['hrrainmat']          #get variable hrrainmat (114x244x4964)
latlist= mat_contents['latlist']              #get varuable latlist (4964)
lonlist= mat_contents['lonlist']


# In[22]:


monthstart=[0,30,61,91,122,153,183,214]
monthend=[30,61,91,122,153,183,214,244]

import numpy as np
april = np.zeros((4964, 114))
may = np.zeros((4964, 114))
june = np.zeros((4964, 114))
july = np.zeros((4964, 114))
aug = np.zeros((4964, 114))
sep = np.zeros((4964, 114))
oct_ber = np.zeros((4964, 114))
nov = np.zeros((4964, 114))


# In[23]:


for i in range(0,4963):    
    april[i,:] = np.mean(hrrainmat[:,monthstart[0]:monthend[0],i], axis=1) 
    may[i,:] = np.mean(hrrainmat[:,monthstart[1]:monthend[1],i], axis=1) 
    june[i,:] = np.mean(hrrainmat[:,monthstart[2]:monthend[2],i], axis=1) 
    july[i,:] = np.mean(hrrainmat[:,monthstart[3]:monthend[3],i], axis=1) 
    aug[i,:] = np.mean(hrrainmat[:,monthstart[4]:monthend[4],i], axis=1) 
    sep[i,:] = np.mean(hrrainmat[:,monthstart[5]:monthend[5],i], axis=1) 
    oct_ber[i,:] = np.mean(hrrainmat[:,monthstart[6]:monthend[6],i], axis=1)
    nov[i,:] = np.mean(hrrainmat[:,monthstart[7]:monthend[7],i], axis=1)


# In[24]:


April = np.zeros((4964, 114,30))
May = np.zeros((4964, 114,31))
June = np.zeros((4964, 114,30))
July = np.zeros((4964, 114,31))
Aug = np.zeros((4964, 114,31))
Sep = np.zeros((4964, 114,30))
Oct = np.zeros((4964, 114,31))
Nov = np.zeros((4964, 114,30))
for i in range(0,4963):    
    April[i,:] = hrrainmat[:,monthstart[0]:monthend[0],i]
    May[i,:] =hrrainmat[:,monthstart[1]:monthend[1],i] 
    June[i,:] = hrrainmat[:,monthstart[2]:monthend[2],i]
    July[i,:] = hrrainmat[:,monthstart[3]:monthend[3],i] 
    Aug[i,:] = hrrainmat[:,monthstart[4]:monthend[4],i]
    Sep[i,:] = hrrainmat[:,monthstart[5]:monthend[5],i] 
    Oct[i,:] = hrrainmat[:,monthstart[6]:monthend[6],i]
    Nov[i,:] = hrrainmat[:,monthstart[7]:monthend[7],i]


# In[25]:


daystart = monthstart[3]
dayend = monthend[6]


# In[26]:


datamat = hrrainmat.transpose(1,0,2)


# In[ ]:





# In[ ]:





# In[27]:


datamatuse = datamat[daystart:dayend,yearstart-1901:yearend-1901,:]


# In[28]:


datamatuse2d = np.squeeze(np.mean(datamatuse , axis = 0))


# In[29]:


datamatuse2d.shape


# In[30]:


precipp = datamatuse2d[:5 , :10]


# In[31]:


precipp.shape


# In[32]:


N = 10


# In[33]:


T = 5


# In[34]:


rainmean = np.zeros((T,1))


# In[35]:


rainmean.shape


# In[36]:


from tqdm import tqdm


# In[37]:


for i in range(0,T) :
    rainmean[i,0] = np.sum(precipp[i,:])


# In[38]:


rainmean[3]


# In[39]:


pmean = np.mean(precipp , axis = 0)


# In[40]:


pmean = np.expand_dims(pmean , axis = 0)


# In[41]:


pmean.shape


# In[42]:


Xp = precipp - np.repeat(pmean ,T , axis = 0)


# In[43]:


Xp.shape


# In[44]:


a = Xp.T


# In[45]:


a.shape


# In[46]:


b = np.dot(a , Xp)


# In[47]:


Sv = b/(T-1)


# In[48]:


Sr = Sv + sigmae*sigmae*np.eye(N)


# In[49]:


varr = np.diag(Sr)


# In[50]:


Er = pmean.T


# In[51]:


options = octave.optimset("Display","iter")


# In[52]:


from ax import RangeParameter, ParameterType


# In[53]:


lst = []

for o in tqdm(range(N)) :
    choice_param = RangeParameter(name=f"location_{o+1}", parameter_type=ParameterType.INT, lower = 0 , upper =1)
    lst += [choice_param]


# In[54]:


lst[8]


# In[55]:


from ax import SumConstraint

con_1 = SumConstraint(parameters= lst, is_upper_bound=False, bound= K)


# In[56]:


con_2 = SumConstraint(parameters= lst, is_upper_bound=True, bound= K)


# In[57]:


from ax import SearchSpace


# In[58]:


search_spc = SearchSpace(parameters=lst, parameter_constraints=[con_1, con_2])


# In[59]:


from ax import *


# In[60]:


def my_evaluation_function(parameterization , weight = None) :
    lst = []
    for g in range(N) :
        lst += [parameterization[f"location_{g+1}"]]
    w = np.expand_dims(np.array(lst) , axis =1)
    beta0 = np.ones((N,1))/np.sum(w)
    ls = np.concatenate((beta0 , w))
    
    betamse = octave.fminunc(octave.getmodelmse,ls , options)
    news = np.concatenate((betamse , w))
    mse = octave.getmodelmse(news)
    return {"objective" : mse}
    
    
    


# In[61]:


exp = SimpleExperiment(
    name="Ws_experiment",
    search_space= search_spc,
    evaluation_function= my_evaluation_function,
    objective_name="objective",
    minimize = True
)


# In[62]:


from ax.modelbridge import get_sobol

sobol = get_sobol(exp.search_space)
exp.new_batch_trial(generator_run=sobol.gen(5))


# In[ ]:





# In[63]:




sobol = Models.SOBOL(search_space=exp.search_space)
generator_run = sobol.gen(5)

for arm in generator_run.arms:
    print(arm)


# In[64]:


exp.new_batch_trial(generator_run=generator_run)


# In[65]:


for arm in exp.trials[0].arms:
    print(arm)


# In[66]:


class MyRunner(Runner):
    def run(self, trial):
        return {"name": str(trial.index)}
    
exp.runner = MyRunner()


# In[67]:


exp.trials[0].run()


# In[ ]:


data = exp.fetch_data()
gpei = Models.BOTORCH(experiment=exp, data=data)
generator_run = gpei.gen(5)
exp.new_batch_trial(generator_run=generator_run)


# In[ ]:





# In[69]:



    
    


# In[70]:



#k = np.array(lj , dtype = 'object')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




