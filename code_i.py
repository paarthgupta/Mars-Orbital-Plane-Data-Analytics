import numpy as np
import pandas as pd
import math
import math
from scipy.optimize import minimize
from scipy.stats.mstats import gmean

#loss :- 0.004097728229404574  
#optimized variables x = 0.9680093969183833   
#optimized variables y = 2.5983627461290855 

df1 = pd.read_csv("./../data/01_data_mars_opposition.csv")
total_data= df1[['Day','Month','Year','ZodiacIndex','Degree','Minute','Second','LatDegree','LatMinute','ZodiacIndexAverageSun','DegreeMean','MinuteMean','SecondMean']]
day=total_data['Day'].values
Month=total_data['Month'].values
Year=total_data['Year'].values
ZodiacIndex=total_data['ZodiacIndex'].values
Degree=total_data['Degree'].values
Minute=total_data['Minute'].values
Second=total_data['Second'].values
ZodiacIndexAverageSun=total_data['ZodiacIndexAverageSun'].values
DegreeMean=total_data['DegreeMean'].values
MinuteMean=total_data['MinuteMean'].values
SecondMean=total_data['SecondMean'].values
alpha=[]
beta=[]
# calculating the alphas and betas
for i in range (12):
  a=ZodiacIndex[i]
  b=Degree[i]
  c=Minute[i]
  d=Second[i]
  alpha.append( a*30 + b + c/60 + d/3600)
  e=ZodiacIndexAverageSun[i]
  f=DegreeMean[i]
  g=MinuteMean[i]
  h=SecondMean[i]
  beta.append( e*30 + f + g/60 + h/3600)
alpha_rad=[]
beta_rad=[]
#converting alphas and betas to radians 
for i in range (12):
  alpha_rad.append(alpha[i] * math.pi / 180.0)
  beta_rad.append(beta[i] * math.pi / 180.0)

var=[1.2,2]

def minimize_distance(var, args):
  r=[]
  alpha_rad = args[0]
  beta_rad = args[1]
  for i in range(12):
    val2 = np.sin(alpha_rad[i] - var[1])
    val1 = var[0] * np.sin(beta_rad[i] - var[1])
    val3 = np.cos(alpha_rad[i] - beta_rad[i])
    radius = np.sqrt((2*val3*val2*val1 + val1*val1 + val2*val2)/(1 - val3*val3))
    r.append(radius)
  ari_mean = np.log(np.mean(r))
  geo_mean = np.log(gmean(r))
  return ari_mean - geo_mean

bnds = ((0, None), (None, None))
out = minimize(minimize_distance, var, args=[alpha_rad,beta_rad],method='L-BFGS-B',bounds=bnds)
var_opt=out['x']
loss=out['fun']
print("loss :- " + str(loss))
print("optimized variables x = " + str(var_opt[0]))
print("optimized variables y = " + str(var_opt[1]))