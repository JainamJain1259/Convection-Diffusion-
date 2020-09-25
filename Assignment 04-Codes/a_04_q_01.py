
#=============Boundary Conditions============#
def apply_bc(phi,num):
    phi[0] = 0
    phi[L] = 1
    for i in range(1,num-1):
        phi[i] = 0.5

#============Thomas Algorithm================#
def Thomas_Algorithm(num, dia, upp, low, rhs): 
    dia1    = np.zeros(num)
    rhs1    = np.zeros(num)
    sol     = np.zeros(num)
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


L        =  1                       # Length of the slab
case     = '2'                      # Case for selecting values of rho,diff and velocity
num_mesh = 11                      # Number of mesh points
del_x    = L/(num_mesh-1.0)         # Mesh size (Delta_x)
del_t    = 0.1                      # Time interval Size(Delta_t)

xmesh   = np.zeros(num_mesh)        # Mesh points
for i in range(0, len(xmesh)):      # Initialising an array of mesh points to locate each node
    xmesh[i] = i * del_x

if case == '1':
    rho  = 1
    diff = 1
    vel    = 1
if case == '2':
    rho  = 1
    diff = 0.01
    vel    = 3

F = rho*vel                          # Strength of convection 
D = diff/del_x                       # Diffusion Strain
P = F/D

print(P)

#-------------------CDS Scheme-------------------#

phi_old_cds    = np.zeros(num_mesh)    # Value of phi at previous timestep in cds scheme
phi_new_cds    = np.zeros(num_mesh)    # Value of phi at current timestep in cds scheme

# Diagonal elements of system matrix
d_cds    = np.zeros(num_mesh)        # main diagonal elements
u_cds    = np.zeros(num_mesh)        # upper diagonal
l_cds    = np.zeros(num_mesh)        # lower diagonal

# RHs of the discretized linear algebraic system
f_cds    = np.zeros(num_mesh)

# apply boundary conditions on the initial field
apply_bc(phi_old_cds,num_mesh)
apply_bc(phi_new_cds,num_mesh)

#Initializing the respective coefficients
a_p_0_cds = rho*del_x/del_t
a_e_cds   = D - (F/2)
a_w_cds   = D + (F/2)
a_p_cds   = a_e_cds + a_w_cds + a_p_0_cds

#Initializing the matrix elements
for i in range(0,num_mesh):
    if i == 0:
       d_cds[i]   = 1
       l_cds[i]   = 0
       u_cds[i]   = 0
       f_cds[i] =   0
    elif i == num_mesh-1:
         d_cds[i]   = 1
         l_cds[i]   = 0
         u_cds[i]   = 0
         f_cds[i] =   1
    else:
         l_cds[i] = -a_w_cds
         u_cds[i] = -a_e_cds
         d_cds[i] =  a_p_cds

iter_cds = 0
converged_cds = False
while converged_cds == False: 
      for i in range(1,num_mesh-1):
          f_cds[i] = a_p_0_cds*phi_old_cds[i]                               #Updating the rhs
      phi_new_cds = Thomas_Algorithm(num_mesh,d_cds,u_cds,l_cds,f_cds)      #Applying thomas algorithm
      norm_cds = abs(max(phi_new_cds-phi_old_cds, key=abs))                 #Calculating the norm

      iter_cds = iter_cds + 1
      #print(iter_cds, norm_cds)
      if norm_cds <0.0000001:
         converged_cds = True
      # Soution not converged -- copy current value into old
      phi_old_cds[:] = phi_new_cds[:]                                       #Updating the solution of previous timestep with current timestep 
print(iter_cds)
#-------------------Upwind Scheme-------------------#

phi_old_u      = np.zeros(num_mesh)         # Value of phi at previous timestep in upwind scheme
phi_new_u      = np.zeros(num_mesh)         # Value of phi at current timestep in upwind scheme

# apply boundary conditions on the initial field
apply_bc(phi_old_u,num_mesh)
apply_bc(phi_new_u,num_mesh)

# Diagonal elements of system matrix
d_u    = np.zeros(num_mesh)        # main diagonal elements
u_u    = np.zeros(num_mesh)        # upper diagonal
l_u    = np.zeros(num_mesh)        # lower diagonal

# RHs of the discretized linear algebraic system
f_u      = np.zeros(num_mesh)

#Initializing the respective coefficients
a_p_0_u = rho*del_x/del_t
a_e_u   = D 
a_w_u   = D + F
a_p_u   = a_e_u + a_w_u + a_p_0_u

#Initializing the matrix elements
for i in range(0,num_mesh):
    if i == 0:
       d_u[i]   = 1
       l_u[i]   = 0
       u_u[i]   = 0
       f_u[i] =   0
    elif i == num_mesh-1:
       d_u[i]   = 1
       l_u[i]   = 0
       u_u[i]   = 0
       f_u[i] =   1
    else:
       l_u[i] = -a_w_u
       u_u[i] = -a_e_u
       d_u[i] =  a_p_u
        
iter_u = 0
converged_u = False
while converged_u == False: 
      for i in range(1,num_mesh-1):
          f_u[i] = a_p_0_u*phi_old_u[i]                       #Updating the rhs
      phi_new_u = Thomas_Algorithm(num_mesh,d_u,u_u,l_u,f_u)  #Applying thomas algorithm
      norm_u = abs(max(phi_new_u-phi_old_u, key=abs))         #Calculating the norm

      iter_u = iter_u + 1
      #print(iter_u, norm_u)
      if norm_u <0.0000001:
         converged_u = True
      # Soution not converged -- copy current value into old
      phi_old_u[:] = phi_new_u[:]                             #Updating the solution of previous timestep with current timestep 
      # plot the converged results
print(iter_u)
#------------------Exact Solution--------------------#   
phi_exact      = np.zeros(num_mesh)    
for i in range(0,num_mesh):
  phi_exact[i] = ((np.exp(rho*vel*xmesh[i]/diff)-1)/(np.exp(rho*vel/diff)-1))

#-----------------Graph of CDS scheme and Exact solution------------------#
plt.plot(xmesh,phi_new_cds,'r-o',marker='v',label='CDS Scheme')
plt.plot(xmesh,phi_exact,'g-o', label ='Exact')
plt.title('CDS Scheme v/s Exact Solution')
plt.xlabel('x')
plt.ylabel('Phi(Φ)')
plt.legend()
plt.show()
#plt.savefig('unsteady-diffusion-implicit.pdf

#-----------------Graph of Upwind scheme and Exact solution------------------#
plt.plot(xmesh,phi_new_u,'b-o',marker ='^',label='Upwind Scheme')
plt.plot(xmesh,phi_exact,'g-o', label = 'Exact')
plt.title('Upwind Scheme v/s Exact Solution')
plt.xlabel('x')
plt.ylabel('Phi(Φ)')
plt.legend()
plt.show()
#plt.savefig('unsteady-diffusion-implicit.pdf''
