

#==============Calculation of sum absolute values of row elements of Error matrix============#
def norm_value(Error,x,y):
    row_elements_a = np.zeros(x)
    for c in range(0,x):
      row_elements_b = np.zeros(y)
      row_elements_c = np.zeros(y)
      for e in range(0,y):
        row_elements_c[e] = Error[c][e]
        row_elements_b[e] = row_elements_b[e-1] + np.absolute(row_elements_c[e])
      row_elements_a[c] = row_elements_b[y-1]
    return(row_elements_a)

#=================Thomas Algorithm=======================#
def Thomas_Algorithm(num, dia, upp, low, rhs):    
    sol     = np.zeros(num)            
    dia1    = np.zeros(num)
    rhs1    = np.zeros(num)
    dia1[0] = dia[0]
    rhs1[0] = rhs[0]
    for i in range(1, len(dia)):
        dia1[i] = dia[i] - low[i]*upp[i-1]/dia1[i-1]
        rhs1[i] = rhs[i] - rhs1[i-1]*low[i]/dia1[i-1]
    sol[num-1] = rhs1[num-1]/dia1[num-1]
    for i in range(len(dia)-2,-1,-1):
        sol[i] = (rhs1[i]-upp[i]*sol[i+1])/dia1[i]
    return(sol)

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys
from scipy.integrate import simps


# Set colour interpolation and colour map.
# You can try set it to 10, or 100 to see the difference
# You can also try: colourMap = plt.cm.coolwarm
colorinterpolation = 50
colourMap = plt.cm.jet

L          = 20                   #Length of the slab
H          = 1                    #Height of the slab 
rho        = 1                    #Density
k          = 1                    #Thermal Conductivity 
c          = 100                  #Specific heat capacity
v          = 0                    #Y direction velocity

num_mesh_x = 601                   #Number of mesh points in x direction
num_mesh_y = 31                    #Number of mesh points in y direction
del_x      = L/(num_mesh_x-1)      #Distance between two consecutive mesh points in x direction
del_y      = H/(num_mesh_y-1)      #Distance between two consecutive mesh points in Y direction

xmesh    = np.zeros(num_mesh_x)    # Mesh point in x direction
ymesh    = np.zeros(num_mesh_y)    # Mesh point in y direction

# Compute the location of mesh points
for i in range(0, len(xmesh)): 
    xmesh[i] = i * del_x
for i in range(0,len(ymesh)):
    ymesh[i] = i*del_y

# Velocity array at each y direction
u = np.zeros(len(ymesh))
for i in range(0,len(ymesh)):
    u[i] = 1.5*(1 -  4*((ymesh[i]-0.5)**2))            #Transformation of axis is performed 

#---------------Initialization of diffusion strain and strength of convection-----------#
F_e = np.zeros((num_mesh_x,num_mesh_y))
F_w = np.zeros((num_mesh_x,num_mesh_y))
F_n = np.zeros((num_mesh_x,num_mesh_y))
F_s = np.zeros((num_mesh_x,num_mesh_y))
D_e = np.zeros((num_mesh_x,num_mesh_y))
D_w = np.zeros((num_mesh_x,num_mesh_y))
D_n = np.zeros((num_mesh_x,num_mesh_y))
D_s = np.zeros((num_mesh_x,num_mesh_y))
P_e = np.zeros((num_mesh_x,num_mesh_y))
P_w = np.zeros((num_mesh_x,num_mesh_y))
P_n = np.zeros((num_mesh_x,num_mesh_y))
P_s = np.zeros((num_mesh_x,num_mesh_y))

#-----------Evaluating the diffusion strain and strength of convection and respective peclet number----------#
for i in range(0,len(xmesh)):
  for j in range(0,len(ymesh)):
      F_e[i,j] = rho*c*u[j]*del_y
      F_w[i,j] = rho*c*u[j]*del_y
      F_n[i,j] = rho*c*v*del_x
      F_s[i,j] = rho*c*v*del_x
      D_e[i,j] = k*del_y/del_x
      D_w[i,j] = k*del_y/del_x
      D_n[i,j] = k*del_x/del_y
      D_s[i,j] = k*del_x/del_y
      P_e[i,j] = F_e[i][j]/D_e[i][j]
      P_w[i,j] = F_w[i][j]/D_w[i][j]
      P_n[i,j] = F_n[i][j]/D_n[i][j]
      P_s[i,j] = F_s[i][j]/D_s[i][j]

#-----------initialization of coefficients---------#
a_e = np.zeros((num_mesh_x,num_mesh_y))
a_w = np.zeros((num_mesh_x,num_mesh_y))
a_n = np.zeros((num_mesh_x,num_mesh_y))
a_s = np.zeros((num_mesh_x,num_mesh_y))
a_p = np.zeros((num_mesh_x,num_mesh_y))

#----------Evaluating respective coefficients using hybrid scheme with respective peclet numbers------------#
for i in range(1,len(xmesh)):
  for j in range(0,len(ymesh)):
        if P_e[i][j] < -2:
            a_e[i,j] = - F_e[i][j]
        elif P_e[i][j] >= -2 and P_e[i][j] <= 2:
            a_e[i,j] = D_e[i][j] - (F_e[i][j]/2)
        elif P_e[i][j] > 2:
            a_e[i,j] = 0

        if P_w[i][j] < -2:
            a_w[i,j] = 0
        elif P_w[i][j] >= -2 and P_w[i][j] <= 2:
            a_w[i,j] = D_w[i][j] + (F_w[i][j]/2)
        elif P_w[i][j] > 2:
            a_w[i,j] = F_w[i][j]

        if P_n[i][j] < -2:
            a_n[i,j] = - F_n[i][j]
        elif P_n[i][j] >= -2 and P_n[i][j] <= 2:
            a_n[i,j] = D_n[i][j] - (F_n[i][j]/2)
        elif P_n[i][j] > 2:
            a_n[i,j] = 0

        if P_s[i][j] < -2:
            a_s[i,j] = 0
        elif P_s[i][j] >= -2 and P_s[i][j] <= 2:
            a_s[i,j] = D_s[i][j] + (F_s[i][j]/2)
        elif P_s[i][j] > 2:
            a_s[i,j] = F_s[i][j]

        a_p[i,j] = a_e[i][j] + a_w[i][j] + a_s[i][j] + a_n[i][j]

T_new          = np.zeros((num_mesh_x,num_mesh_y))
T_old          = np.zeros((num_mesh_x,num_mesh_y))
T_bulk         = np.zeros(num_mesh_x)
T_1            = np.zeros(num_mesh_x)
Error          = np.zeros((num_mesh_x,num_mesh_y))

#-------Boudary condition-------#
T_above = 100.0
T_below = 100.0
T_left = 50.0 

for i in range(0,len(xmesh)):
  for j in range(0,len(ymesh)):
    if j == 0:
      T_new[i,j] = T_below
      T_old[i,j] = T_below
    elif j == len(ymesh)-1 :
      T_new[i,j] = T_above
      T_old[i,j] = T_above
    elif i == 0:
      T_new[i,j] = T_left
      T_old[i,j] = T_left
    else:  
      T_old[i,j] = 50 
      T_new[i,j] = 50    
                    
#----------Matrix Elements-----------#
low = np.zeros(num_mesh_x)
upp = np.zeros(num_mesh_x)
dia = np.zeros(num_mesh_x)
rhs = np.zeros(num_mesh_x)
solution = np.zeros(num_mesh_x)
row_elements = np.zeros(num_mesh_x)

iter = 0
while True:
       for j in range(1,len(ymesh)-1):
        for i in range(0,len(xmesh)):
          if i==0:                                     #matrix elements at x=0
            dia[i] = 1
            low[i]= 0
            upp[i]= 0
            rhs[i]= 50
          elif i==num_mesh_x-1:                        #matrix elements at x=len(xmesh)-1
            dia[i] = 1
            low[i] = -1
            upp[i] = 0
            rhs[i] = 0
          else:                                        #matrix elements else where
            low[i] = -a_w[i][j]
            upp[i] = -a_e[i][j]
            dia[i] =  a_p[i][j]
            rhs[i] = a_s[i][j]*T_new[i][j-1] + a_n[i][j]*T_new[i][j+1] 
        solution = Thomas_Algorithm(len(xmesh),dia,upp,low,rhs)
        for b in range(0,len(xmesh)):
            T_new[b,j] = solution[b]
               
       Error[:][:] = T_new[:][:] - T_old[:][:]
       row_elements = norm_value(Error,len(xmesh),len(ymesh))
       norm = np.max(row_elements)
       #print(iter,norm)
       if(norm<0.000001):
         break
       iter = iter + 1
       T_old[:][:] = T_new[:][:]
print("Number of iterations:",iter)

# Configure the contour
plt.title("Contour of Temperature")
plt.contourf(xmesh, ymesh, np.transpose(T_new), colorinterpolation, cmap=colourMap)
# Set Colorbar
plt.colorbar()
# Show the result in the plot window
plt.show()

T_bulk = np.zeros(len(xmesh))             #Bulk mean temperature at each axial location

for i in range(0,len(xmesh)):
    T_1 = np.zeros(len(ymesh))
    for j in range(0,len(ymesh)): 
        T_1[j] = T_new[i][j]*u[j]
    T_bulk[i] = simps(T_1,ymesh)
#print(T_bulk)

h = np.zeros(num_mesh_x)                  #Heat transfer coefficient at each axial location
Nu = np.zeros(num_mesh_x)                 #nusselt number at each axial location

for i in range(0,len(xmesh)):
    h[i] = k*(T_below-T_new[i][1])/((T_below-T_bulk[i])*del_y)
    Nu[i] = h[i]*2*H/k

#print(Nu)

plt.plot(xmesh,Nu,'r-o')
plt.xlabel('Mesh Points')
plt.ylabel('Nusslet Number')
plt.show()
