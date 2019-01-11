#!/usr/bin/env python
# coding: utf-8

# In[1]:


def f1(x):
    return np.sin(np.sqrt(x[0] ** 2 + x[1] ** 2)) * (np.sqrt(x[0] ** 2 + x[1] ** 2))


# In[2]:


import numpy as np
x = np.linspace(-20, 20, num=1000)
y = np.linspace(-20, 20, num=1000)
#print(x)

func=f1
f_list=[]
for x_val in x:
    for y_val in y:
        x_input = np.array([x_val,y_val])
        f_list.append(func(x_input))
print("Min / Max: ", min(f_list),max(f_list))


# In[3]:


from scipy import optimize
np.random.seed(555)   # Seeded to allow replication.

# In[5]:


bounds = [(-20, 20), (-20, 20)]
de = optimize.differential_evolution(f1,bounds)
de.x, de.fun


# In[6]:


da = optimize.dual_annealing(f1,bounds=bounds)
da.x, da.fun


# In[7]:


import LengthToGamma as lg
'''
lengths_list = [[143, 313, 328, 135, 132, 167, 165],
                [181, 152, 307, 100, 259, 100, 199],
                [177, 145, 113, 100, 259, 100, 100],
                [174, 145, 316, 100, 267, 100, 389],
                [160, 340, 299, 95, 281, 299, 280],
                250 * np.ones(20),
                260 * np.ones(20),
                270 * np.ones(20),
                280 * np.ones(20),
                290 * np.ones(20)]

for lengths in lengths_list:
    print(lg.LengthToGamma(lengths))
'''
def optGamma(lengths):
    return -lg.LengthToGamma(lengths)

f1 = optGamma

# In[15]:

bounds = [(10, 56),(10, 56),(10, 56),(10, 56),(10, 56),(10, 56),(10, 56),(10, 56),(10, 56),(10, 56)]
#bounds = [(40,400),(40,400),(40,400)]
da = optimize.dual_annealing(f1,bounds=bounds,maxiter=1000)
da.x, -da.fun

#bounds = [(40,400),(40,400),(40,400),(40,400),(40,400),(40,400),(40,400),(40,400),(40,400),(40,400),(40,400),(40,400),(40,400)]
#da = optimize.dual_annealing(f1,bounds=bounds,maxiter=5000)
#da.x, -da.fun

# In[13]:


da.message


# In[14]:


vals = da.x.astype(int)
nm_vals = [10 * val for val in vals] #Conversion from 10*nm to nm (see note in LengthToGamma function in LengthToGamma.py)
print("Lengths Vector: ", nm_vals)
print("Corresponding Gamma: ", -f1(list(vals)))


# In[ ]:
