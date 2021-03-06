

#==============Calculation of sum absolute values of row elements of Error matrix============#
def norm_value(Error,x,y):
    for c in range(0,x):
      row_elements = np.zeros(y)
      row_elements_1 = np.zeros(y)
      for e in range(0,y):
        row_elements_1[e] = Error[c][e]
        row_elements[e] = row_elements[e-1] + np.absolute(row_elements_1[e])
      row_elements_2[c] = row_elements[y-1]
    return(row_elements_2)

#=================Thomas Algorithm=======================#
def Thomas_Algorithm(num, dia, upp, low, rhs):                
    dia1    = np.zeros(num)
    rhs1    = np.zeros(num)
    dia1[0] = dia[0]
    rhs1[0] = rhs[0]
    for i in range(1, len(dia)):
        dia1[i] = dia[i] - low[i]*upp[i-1]/dia1[i-1]
        rhs1[i] = rhs[i] - rhs1[i-1]*low[i]/dia1[i-1]
    sol[num-1] = rhs[num-1]/dia[num-1]
    for i in range(len(dia)-2,-1,-1):
        sol[i] = (rhs1[i]-upp[i]*sol[i+1])/dia1[i]
    return(sol)

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

# Set colour interpolation and colour map.
# You can try set it to 10, or 100 to see the difference
# You can also try: colourMap = plt.cm.coolwarm
colorinterpolation = 50
colourMap = plt.cm.jet

L          = 20                    #Length of the slab
H          = 20                    #Height of the slab
scheme     = 'CDS'                 #Scheme to be used for discretizations
Method     = 'Line by Line'        #Method to be used to solve the equation
num_mesh_x = 21                    #Number of mesh points in x direction
num_mesh_y = 21                    #Number of mesh points in y direction
del_x      = L/(num_mesh_x-1)      #Distance between two consecutive mesh points in x direction
del_y      = H/(num_mesh_y-1)      #Distance between two consecutive mesh points in Y direction

xmesh    = np.zeros(num_mesh_x)        # Mesh point in x direction
ymesh    = np.zeros(num_mesh_y)        # Mesh point in y direction

# Compute the location of mesh points
for i in range(0, len(xmesh)):
    xmesh[i] = i * del_x
for i in range(0, len(ymesh)):
    ymesh[i] = i * del_y

#---------Given Values in the question--------#
rho  = 1
diff = 1
u    = 1
v    = 4
a    = 10
b    = 2

#---------Defining the variables----------#
phi_new = np.zeros((num_mesh_x,num_mesh_y))
phi_old = np.zeros((num_mesh_x,num_mesh_y))
Error = np.zeros((len(xmesh), len(ymesh)))
row_elements_2 = np.zeros(num_mesh_y)

#----------Boundary Conditions-------#
phi_above = 0.0
phi_below = 100.0
phi_right = 0.0
phi_left  = 100.0
phi_guess = 0.0

for i in range(0,len(xmesh)):
  for j in range(0,len(ymesh)):
    if j == 0:
      phi_new[i,j] = phi_below
      phi_old[i,j] = phi_below
    elif j == len(ymesh)-1 :
      phi_new[i,j] = phi_above
      phi_old[i,j] = phi_above
    elif i == 0:
      phi_new[i,j] = phi_left
      phi_old[i,j] = phi_left
    elif i == len(xmesh)-1:
      phi_new[i,j] = phi_right
      phi_old[i,j] = phi_right
    else:  
      phi_old[i,j] = phi_guess

#---------Calculation of diffusive strain and Strength of convection---------#
D_e = diff*del_y/del_x
F_e = rho*u*del_y

D_w = diff*del_y/del_x
F_w = rho*u*del_y

D_n = diff*del_x/del_y
F_n = rho*v*del_x

D_s = diff*del_x/del_y
F_s = rho*v*del_x

#-----------Defining the respective Peclet Numbers-------------#
P_e = F_e/D_e
P_w = F_w/D_w
P_s = F_s/D_s
P_n = F_n/D_n

#----------Assigning the values to the coefficients according to the scheme chosen------------#
if scheme == 'CDS':
    a_w  = D_w + (F_w/2)
    a_e  = D_e - (F_e/2)
    a_s  = D_s + (F_s/2)
    a_n  = D_n - (F_n/2)
    a_p  = a_w + a_e + a_s + a_n + b*del_x*del_y
    l    = a*del_y*del_x

elif scheme == 'Upwind':
    a_w  = D_w + F_w
    a_e  = D_e
    a_s  = D_s + F_s
    a_n  = D_n
    a_p  = a_w + a_e + a_s + a_n + b*del_x*del_y
    l    = a*del_y*del_x

elif scheme == 'Hybrid':

    if P_e < -2:
        a_e = - F_e
    elif P_e >= -2 and P_e <= 2:
        a_e = D_e - (F_e/2)
    elif P_e >= 2:
        a_e = 0

    if P_w < -2:
        a_w = 0
    elif P_w >= -2 and P_w <= 2:
        a_w = D_w + (F_w/2)
    elif P_w >= 2:
        a_w = F_w

    if P_s < -2:
        a_n = - F_n
    elif P_n >= -2 and P_n <= 2:
        a_n = D_n - (F_n/2)
    elif P_n >= 2:
        a_n = 0

    if P_s < -2:
        a_s = 0
    elif P_s >= -2 and P_s <= 2:
        a_s = D_s + (F_s/2)
    elif P_s >= 2:
        a_s = F_s
    
    a_p  = a_w + a_e + a_s + a_n + b*del_x*del_y
    l    = a*del_y*del_x

#--------------Point by Point Method------------#
if Method == 'Point by Point':
            converged = False
            iter = 0

            while converged == False:
              iter = iter + 1

              for i in range(1,len(xmesh)-1):
                for j in range(1,len(ymesh)-1):
                  phi_new[i,j] = (a_w*phi_new[i-1][j] + a_e*phi_new[i+1][j] + a_s*phi_new[i][j-1] + a_n*phi_new[i][j+1] + l)/a_p
                
              Error[:][:] = phi_new[:][:]- phi_old[:][:]
              row_elements_2 = norm_value(Error,len(xmesh),len(ymesh))

              norm = np.max(row_elements_2) 

              if(norm<0.00000001):
                converged = True

              phi_old[:][:] = phi_new[:][:] 
            print("Number of iterations:",iter)
            # Configure the contour
            plt.title(str(Method) + " Method and " + str(scheme) + " Scheme")
            plt.contourf(xmesh, ymesh, np.transpose(phi_new), colorinterpolation, cmap=colourMap)

            # Set Colorbar
            plt.colorbar()

            # Show the result in the plot window
            plt.show()

#----------Matrix Elements-----------#
low = np.zeros(num_mesh_y)
upp = np.zeros(num_mesh_y)
dia = np.zeros(num_mesh_y)
rhs = np.zeros(num_mesh_y)
sol = np.zeros(num_mesh_y)

#=========Initializing the matrix elements===========#

for j in range(0,num_mesh_y):
  if j == 0:
    dia[j]   = 1
    low[j]   = 0
    upp[j]   = 0
    rhs[j]   = 100
  elif j == num_mesh_y-1:
    dia[j]   = 1
    low[j]   = 0
    upp[j]   = 0
    rhs[j]   = 0
  else:
    low[j] = -a_s
    upp[j] = -a_n
    dia[j] =  a_p

iter = 0
#-----------Line by Line Method------------#
if Method == 'Line by Line':

            while True:
                  for i in range(1,len(xmesh)-1):
                      for j in range(1,len(ymesh)-1):
                          rhs[j] = a_w*phi_new[i-1][j] + a_e*phi_new[i+1][j] + l 
                      sol = Thomas_Algorithm(num_mesh_y,dia,upp,low,rhs)
                      for b in range(0,num_mesh_y):
                          phi_new[i,b] = sol[b]
                  Error[:][:] = phi_new[:][:] - phi_old[:][:]
                  row_elements_2 = norm_value(Error,len(xmesh),len(ymesh))
                  norm = np.max(row_elements_2)

                  if(norm<0.000001):
                    break
                  iter = iter + 1
                  phi_old[:][:] = phi_new[:][:]
            print("Number of iterations:",iter) 
            # Configure the contour
            plt.title(str(Method) + " Method and " + str(scheme) + " Scheme")
            plt.contourf(xmesh, ymesh, np.transpose(phi_new), colorinterpolation, cmap=colourMap)

            # Set Colorbar
            plt.colorbar()

            # Show the result in the plot window
            plt.show()
