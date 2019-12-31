import numpy as np
import pandas as pd
import math
import math
from scipy.optimize import minimize
from scipy.stats.mstats import gmean

#Five different projections of Mars location on the ecliptic plane (in radians) are :
#[2.604342383281199, 5.761762953132912, 0.7746792940121772, 3.0508413387651125, 2.759205698738966]
#loss :- 0.07120993396639322                          
#Radius in AU = 1.5773209107563353

df1 = pd.read_csv("./../data/01_data_mars_opposition.csv")
df2 = pd.read_csv("./../data/01_data_mars_triangulation.csv")
total_data= df1[['Day','Month','Year','ZodiacIndex','Degree','Minute','Second','LatDegree','LatMinute','ZodiacIndexAverageSun','DegreeMean','MinuteMean','SecondMean']]
day=total_data['Day'].values
Month=total_data['Month'].values
Year=total_data['Year'].values
ZodiacIndex=total_data['ZodiacIndex'].values
Degree=total_data['Degree'].values
Minute=total_data['Minute'].values
Second=total_data['Second'].values
Latdegree=total_data['LatDegree'].values
LatMinute=total_data['LatMinute'].values
ZodiacIndexAverageSun=total_data['ZodiacIndexAverageSun'].values
DegreeMean=total_data['DegreeMean'].values
MinuteMean=total_data['MinuteMean'].values
SecondMean=total_data['SecondMean'].values

total_data2= df2[['PairIndex','Day','Month','Year','DegreeEarthLocationHelioCentric','MinuteEarthLocationHelioCentric','DegreeMarsLocationGeoCentric','MinuteMarsLocationGeoCentric']]
Year2=total_data2['PairIndex'].values
day2=total_data2['Day'].values
Month2=total_data2['Month'].values
DegreeEarthLocationHelioCentric2=total_data2['DegreeEarthLocationHelioCentric'].values
MinuteEarthLocationHelioCentric2=total_data2['MinuteEarthLocationHelioCentric'].values
DegreeMarsLocationGeoCentric2=total_data2['DegreeMarsLocationGeoCentric'].values
MinuteMarsLocationGeoCentric2=total_data2['MinuteMarsLocationGeoCentric'].values

alpha=[]
beta=[]
# calculating the alphas and betas
for i in range (10):
  b=DegreeEarthLocationHelioCentric2[i]
  c=MinuteEarthLocationHelioCentric2[i]
  alpha.append( b + c/60)
  f=DegreeMarsLocationGeoCentric2[i]
  g=MinuteMarsLocationGeoCentric2[i]
  beta.append( f + g/60)
alpha_rad=[]
beta_rad=[]
#converting alphas and betas to radians 
for i in range (10):
  alpha_rad.append(alpha[i] * math.pi / 180.0)
  beta_rad.append(beta[i] * math.pi / 180.0)

r=[]
mars_pos=[]
x_coordinates_list=[]
y_coordinates_list=[]
for i in range (5):
  theta1=alpha_rad[2*i]
  theta2=alpha_rad[2*i + 1]
  alpha1=beta_rad[2*i]
  alpha2=beta_rad[2*i + 1]
  tan_alpha1=np.tan(alpha1)
  tan_alpha2=np.tan(alpha2)
  den=tan_alpha1 - tan_alpha2
  x= ((-1 / den) * (np.sin(theta1)-(tan_alpha1 * np.cos(theta1)) )) +(1 / den) * (np.sin(theta2)-(tan_alpha2 * np.cos(theta2)) )
  y= ((-tan_alpha2  / den) * (np.sin(theta1)-(tan_alpha1 * np.cos(theta1)) )) +(tan_alpha1 / den) * (np.sin(theta2)-(tan_alpha2 * np.cos(theta2)) )
  r.append(np.sqrt(x*x + y*y))  
  # for different cases of x and y we calculate the angle
  if x > 0 and y > 0:
      angle = np.arctan(y / x)
  elif x > 0 > y:
      angle = 2 * math.pi - np.arctan(np.abs(y) / x)
  elif y > 0 > x:
      angle = math.pi - np.arctan(y / np.abs(x))
  else:
      angle = math.pi - np.arctan(np.abs(y) / np.abs(x))
	  
  x_coordinates_list.append(x)
  y_coordinates_list.append(y)
  mars_pos.append(angle)
  
print("Five different projections of Mars location on the ecliptic plane (in radians) are :") 
print(mars_pos)
print("x: " + str(x_coordinates_list))
print("y: " + str(y_coordinates_list))

def minimize_distance(var, args):
  dis=[]
  r = var[0]
  d = args[0]
  for i in range(5):
      z = math.pow(r - d[i], 2)
      dis.append(z)
  return np.sum(dis)
 
  
  
var=[1.5]
bnds = [(1.2, None)]
out = minimize(minimize_distance, var, args=[r],method='L-BFGS-B',bounds=bnds)
var_opt=out['x']
loss=out['fun']
print("loss :- " + str(loss))
print("Radius in AU = " + str(var_opt[0]))