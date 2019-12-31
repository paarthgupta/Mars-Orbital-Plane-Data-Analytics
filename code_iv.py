import numpy as np
import pandas as pd
import math
import math
from scipy.optimize import minimize
from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt 
import matplotlib.patches as pt

#coordinates of mars in 3 dimentions
#x: [-1.4529736727603795, 1.195672782788594, 1.0738853142069968, -1.632304590013056, -1.553767331486135]    
#y: [0.8655335301531041, -0.6868566346181093, 1.0511069275483509, -0.14854179871578346, 0.6248989852957588] 
#z: [0.053483392412438116, -0.04345068547829168, -0.0021522535928766547, 0.03530594693180208, 0.050546831728016724] 
#loss :- 0.07132672368470987                 
#Radius in AU = 1.5778557630066357
#loss :- 0.019689864425727954
#x coordinate of focus 2 :- -0.24394173355009974  
#y coordinate of focus 2 :- 0.20144851243890732                            
#length of major axis :- 3.0706447265643626                                                        
#loss function used is minimizing the ((sum of (sum of distances from a given point and the 2 foci of that ellipse)-length of major axis)) 

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
z_coordinates_list=[]
inc=1.8538035741165297

for i in range (5):
  theta1=alpha_rad[2*i]
  theta2=alpha_rad[2*i + 1]
  alpha1=beta_rad[2*i]
  alpha2=beta_rad[2*i + 1]
  tan_alpha1=np.tan(alpha1)
  tan_alpha2=np.tan(alpha2)
  den=tan_alpha1 - tan_alpha2
  x= ((-1 / den) * (np.sin(theta1)-(tan_alpha1 * np.cos(theta1)))) +(1 / den) * (np.sin(theta2)-(tan_alpha2 * np.cos(theta2)) )
  y= ((-tan_alpha2  / den) * (np.sin(theta1)-(tan_alpha1 * np.cos(theta1)) )) +(tan_alpha1 / den) * (np.sin(theta2)-(tan_alpha2 * np.cos(theta2)) )
  x_coordinates_list.append(x)
  y_coordinates_list.append(y)
  r.append(np.sqrt(x*x + y*y))  
    
  if x > 0 and y > 0:
      angle = np.arctan(y / x)
  elif x > 0 > y:
      angle = 2 * math.pi - np.arctan(np.abs(y) / x)
  elif y > 0 > x:
      angle = math.pi - np.arctan(y / np.abs(x))
  else:
      angle = math.pi - np.arctan(np.abs(y) / np.abs(x))

  mars_pos.append(angle)
a3=0.6470968384614812
b3=-0.6050737005314765
c3=27.371604156028972
#print(mars_pos)
for i in range(5):
  z_coordinates_list.append( ((-a3*x_coordinates_list[i])- b3*y_coordinates_list[i])/c3)
print("coordinates of mars in 3 dimentions")
print("x: " + str(x_coordinates_list))
print("y: " + str(y_coordinates_list))
print("z: " + str(z_coordinates_list))

distance=[]
for i in range(5):
  distance.append(np.sqrt(x_coordinates_list[i] * x_coordinates_list[i] + y_coordinates_list[i] * y_coordinates_list[i] + z_coordinates_list[i] * z_coordinates_list[i]))


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
out = minimize(minimize_distance, var, args=[distance],method='L-BFGS-B',bounds=bnds)
var_opt=out['x']
loss=out['fun']
print("loss :- " + str(loss))
print("Radius in AU = " + str(var_opt[0]))

c=plt.Circle((0,0),1.5773209107563353,fill=False)
fig,ax=plt.subplots(figsize=(10,10))
ax.add_artist(c)
plt.scatter(x_coordinates_list,y_coordinates_list)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.show()
#plt.savefig("2_iv_fig1")

def minimize_distance(var, args):
    coordinate_x=args[0]
    coordinate_y=args[1]

    dis=[]
    for i in range(5):
      
      zz= abs(np.sqrt(coordinate_x[i]*coordinate_x[i] + coordinate_y[i]*coordinate_y[i]) + np.sqrt(np.square(coordinate_x[i]-var[0])+np.square(coordinate_y[i]-var[1])) - var[2])
         
      dis.append(zz)
    #print(dis)
    #print(np.sum(dis)) 
    #print(np.sum(dis))
    return np.sum(dis)  
var=[1,1,1]
#bnds = [(1.2, None)]
out = minimize(minimize_distance, var, args=[x_coordinates_list,y_coordinates_list],method='L-BFGS-B')
var_opt=out['x']
#inc=np.arccos(var_opt[2]/np.linalg.norm(var_opt))*(180/np.pi)
loss=out['fun']
print("loss :- " + str(loss))
print("x coordinate of focus 2 :- " + str(var_opt[0]))
print("y coordinate of focus 2 :- " + str(var_opt[1]))
print("length of major axis :- " + str(var_opt[2]))
print("loss function used is minimizing the ((sum of (sum of distances from a given point and the 2 foci of that ellipse)-length of major axis))")

x=-0.24394173/2
y=0.20144851/2
c=np.sqrt(np.square(-0.24394173) + np.square(0.20144851)) / 2
b=np.sqrt(np.square(3.07064473 / 2)- np.square(c) )
d=pt.Ellipse( (x,y),3.07064473,2 *b   ,angle=0,fill=False)
fig,ax=plt.subplots(figsize=(10,10))
ax.add_artist(d)
plt.scatter(x_coordinates_list,y_coordinates_list)

xx=[0,x]
yy=[0,y]
plt.scatter(xx,yy)

plt.xlim([-2,2])
plt.ylim([-2,2])
plt.show()
#plt.savefig("2_iv_fig2")

y4=[1.4394458945438964,1.5014896753783606,0.9747386019545708,0.11722917382644234,-0.8865956950779831,-1.566917999583486,-0.47785030102239034,1.1580390999650196,1.5681475650401604,1.2259167541156522,0.4732490965051435,-0.5012114194578124]

x4=[0.6265943464442608,-0.4569031621546369,-1.2299365359470285,-1.5651874087757587,-1.2956486880860487,-0.08989061772290922,1.4943008401709512,1.0601149388047009,-0.06755277403225604,-0.9798047440522887,-1.4963997424804596,-1.4876493315653687]

c=plt.Circle((0,0),1.5773209107563353,fill=False)
fig,ax=plt.subplots(figsize=(10,10))
ax.add_artist(c)
plt.scatter(x4,y4)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.show()
#plt.savefig("2_iv_fig3")