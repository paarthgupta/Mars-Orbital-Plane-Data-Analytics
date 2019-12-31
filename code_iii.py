import numpy as np
import pandas as pd
import math
import math
from scipy.optimize import minimize
from scipy.stats.mstats import gmean

#The corresponding heliocentric latitudes of Mars are
#[0.0106495055508308, 0.026230164044516092, 0.029012019659864406, 0.023557766834431733, 0.0076667418302551775, -0.025588578049221964, -0.03866560323233603, 0.0008517536638064522, 0.022703068778338335, 0.02879793178276269, 0.026657965022637908, 0.015552578445959982]
#The corresponding heliocentric longitude of Mars are
#[1.160231860947286, 1.866193302694937, 2.4714347021920813, 3.0668343839627026, 3.741695028067177, 4.655084003277542, 5.973680253159258, 0.829516208378416, 1.6138477816774235, 2.2450751944820393, 2.8352873698647882, 3.466563264037515]
#loss :- 0.0003649683032666839
#Inclination = 1.8538035741165297
#Parameters of the plane are
#a:0.6470968384614812
#b:-0.6050737005314765
#c:27.371604156028972

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

calculated_radius=1.5773209107563353
#calculated_radius=1.57
alpha=[]
alpha_rad_lat=[]
lat_mar_heli=[]
coordinate_x = []
coordinate_y = []
coordinate_z = []
#beta=[]
# calculating the alphas and betas
for i in range (12):
  a=ZodiacIndex[i]
  b=Degree[i]
  c=Minute[i]
  d=Second[i]
  alpha.append( a*30 + b + c/60 + d/3600)
  #e=ZodiacIndexAverageSun[i]
  #f=DegreeMean[i]
  #g=MinuteMean[i]
  #h=SecondMean[i]
  #beta.append( e*30 + f + g/60 + h/3600)
alpha_rad=[]
#beta_rad=[]
#converting alphas and betas to radians 
for i in range (12):
  alpha_rad.append(alpha[i] * math.pi / 180.0)
  #beta_rad.append(beta[i] * math.pi / 180.0)
for i in range (12):
  alpha_rad_lat.append((math.pi / 180.0)* (Latdegree[i] + LatMinute[i] /60))
for i in range (12):
  lat_mar_heli.append(np.arctan(np.tan(alpha_rad_lat[i])*((calculated_radius-1)/calculated_radius)))
# let the plane equation be A1x + A2y + A3z + d =0
coordinate_z=calculated_radius*np.sin(lat_mar_heli)
coordinate_x=calculated_radius*np.cos(lat_mar_heli)*np.cos(alpha_rad)
coordinate_y=calculated_radius*np.cos(lat_mar_heli)*np.sin(alpha_rad)
#print(coordinate_x)
#print(coordinate_y)
#print(coordinate_z)
print("The corresponding heliocentric latitudes of Mars are")
print(lat_mar_heli)
print("The corresponding heliocentric longitude of Mars are")
print(alpha_rad)

def minimize_distance(var, args):
    coordinate_x=args[0]
    coordinate_y=args[1]
    coordinate_z=args[2]
    dis=[]
    for i in range(12):
      #print(i)
      #print(coordinate_x)
      a=abs((var[0] * coordinate_x[i] + var[1] * coordinate_y[i] + var[2] * coordinate_z[i]  ))
      #print(a)
      b=(math.sqrt(var[0] * var[0] + var[1] * var[1] + var[2] * var[2]))
      #print(b)
      
      dis.append(np.square(a / b))
    #print(dis)
    #print(np.sum(dis))  
    return np.sum(dis)  
var=[10,10,10]
#bnds = [(1.2, None)]
out = minimize(minimize_distance, var, args=[coordinate_x,coordinate_y,coordinate_z],method='L-BFGS-B')
var_opt=out['x']
aa=var_opt[0]
bb=var_opt[1]
cc=var_opt[2]
#dd=var_opt[3]
inc=np.arccos(var_opt[2]/np.linalg.norm(var_opt))*(180/np.pi)
loss=out['fun']
print("loss :- " + str(loss))
print("Inclination in degrees= " + str(inc))
#print("Radius in AU = " + str(var_opt[0]))
print("Parameters of the plane are")
print("a:" + str(aa))
print("b:" + str(bb))
print("c:" + str(cc))
#print("d:" + str(dd))





