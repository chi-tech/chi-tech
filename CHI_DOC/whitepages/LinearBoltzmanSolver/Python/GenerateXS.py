#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 20:03:14 2018

@author: janv4
"""
import numpy as np
from pyne import data as pyne_data
from pyne import ace  as pyne_ace
import matplotlib.pyplot as plt
import math
import Scattering

#============================ Define Isotope and ACE file
A=pyne_data.atomic_mass('C12')
Acefile='6000.710nc'
AcePath='/Users/janv4/Desktop/Projects/' + \
        'MCNPDATA/MCNP_DATA/DVD-3/' + \
        'MORE_DATA_620-B/xdata/endf71x/C/'
OutputFileName='SXS_6000.txt'
OutputPath='/Users/janv4/Desktop/Personal/' + \
           'compphysicsnotes/NeutronTransport/SXS/Groups_10A/'
           
OF = open(OutputPath+OutputFileName,'w')


           
#============================ Define flux spectrum to use as weighting
FluxWeightOption='Unity'
        
#============================ Define Groups
G = 10                              #Number of energy groups
Eg = np.array([10.0,2.0,1.5, \
               1.1,0.9,0.2, \
               0.01,100e-6,10e-6,1e-6,1e-10])
int_res = 1000                      #Integration resolution
print("Energy group boundaries")
print(Eg)

OF.write('# Number of groups='+str(G)+'\n')
for g in range(0,G+1):
    OF.write(str('%.4e \n' %Eg[g]))

#============================ Define max scattering order
Lmax = 9

#============================ Plot properties
plot_xmin=1e-10
plot_xmax=200
plot_ymin=1e-6
plot_ymax=1e6

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# DONT CHANGE ANYTHING BELOW THIS
        
#============================ Extracting ACE data
print("A=%f" %A)
pyne_libfile = pyne_ace.Library(AcePath+Acefile)
pyne_libfile.read()

pyne_XSTable = pyne_libfile.tables[Acefile]

print('The following reactions are available:')
print('<ACE Reaction: MT=1 >')
for MT in pyne_XSTable:
    print(MT)
    
pyne_E = pyne_XSTable.energy
pyne_sigma_t = pyne_XSTable.sigma_t
pyne_sigma_s = pyne_XSTable.reactions[2].sigma
pyne_sigma_a = pyne_XSTable.reactions[102].sigma

out_sigma_t  = np.zeros((G))
out_sigma_s  = np.zeros((G))
out_sigma_a  = np.zeros((G))

#============================ Compute multi-group cross-sections
for g in range(0,G):
    dE = (Eg[g]-Eg[g+1])/int_res
    #EfineArr = np.linspace(Eg[g],Eg[g+1],int_res)
    #EfineArr = np.logspace(math.log10(Eg[g+1]),math.log10(Eg[g]),int_res)
    sum_flux=0
    sum_xs_t = 0
    sum_xs_s = 0
    sum_xs_a = 0
    for gi in range(0,int_res):
        #Efine = EfineArr[gi]
        Efine = Eg[g+1] + 0.5*dE + dE*gi
        
        if (FluxWeightOption=='Unity'):
            fluxperunitE = 1.0
            
        sum_xs_t = sum_xs_t + \
            np.interp(Efine,pyne_E,pyne_sigma_t)* \
            fluxperunitE*1.0*dE
        sum_xs_s = sum_xs_s + \
            np.interp(Efine,pyne_E,pyne_sigma_s)* \
            fluxperunitE*1.0*dE
        sum_xs_a = sum_xs_a + \
            np.interp(Efine,pyne_E,pyne_sigma_a)* \
            fluxperunitE*1.0*dE
        sum_flux = sum_flux + 1.0*dE
    
    out_sigma_t[g] = sum_xs_t/sum_flux
    out_sigma_s[g] = sum_xs_s/sum_flux
    out_sigma_a[g] = sum_xs_a/sum_flux

#============================ Print Sigma_t,s and a to file
OF.write('# MT1, total ==================================\n')
k=0
for g in range(0,G):
   k=k+1
   if (out_sigma_t[g]<0.001):
       OF.write('%14.4e ' %out_sigma_t[g])
   else:
       OF.write('%14.4f ' %out_sigma_t[g])
   if (k==5):
       k=0
       OF.write('\n')

OF.write('# MT2, elastic scattering ===============================\n')
k=0
for g in range(0,G):
   k=k+1
   if (out_sigma_s[g]<0.001):
       OF.write('%14.4e ' %out_sigma_s[g])
   else:
       OF.write('%14.4f ' %out_sigma_s[g])
   if (k==5):
       k=0
       OF.write('\n')       

OF.write('# MT102, z,gamma =======================================\n')
k=0
for g in range(0,G):
   k=k+1
   if (out_sigma_a[g]<0.001):
       OF.write('%14.4e ' %out_sigma_a[g])
   else:
       OF.write('%14.4f ' %out_sigma_a[g])
   if (k==5):
       k=0
       OF.write('\n')  
       
#============================ Make stair-step array
StepArrSize = ((G+1)-2)*2 +2
StepE    = np.zeros((  StepArrSize      ))
StepXS_t = np.zeros((  StepArrSize      ))
StepXS_s = np.zeros((  StepArrSize      ))
StepXS_a = np.zeros((  StepArrSize      ))
b  = 0
gi = 0
k  = 0
StepE[k] = Eg[b]
StepXS_t[k] = out_sigma_t[gi]
StepXS_s[k] = out_sigma_s[gi]
StepXS_a[k] = out_sigma_a[gi]

flipflop = 0
while (k<(StepArrSize-1)):
    k=k+1
    if (flipflop == 0):
        flipflop = 1
        b=b+1
    else:
        flipflop = 0
        gi = gi +1
    StepE[k]    = Eg[b]
    StepXS_t[k] = out_sigma_t[gi]
    StepXS_s[k] = out_sigma_s[gi]
    StepXS_a[k] = out_sigma_a[gi]
        
#============================ Compute Transfer matrices
L=Lmax
for ell in range(0,L+1):
    OF.write('# Transfer matrix Moment %d ======================\n' %ell)
    for g in range(0,G):
        OF.write('# Scattering to group %d\n' %g)
        k=0
        for gprime in range(0,G):
            k=k+1
            Kl = out_sigma_s[gprime]*Scattering.KernelL(ell,gprime,g,Eg,A)
            print("Moment %d, %d->%d, %f)" %(ell,gprime,g,Kl))
            if (Kl<0.001):
                OF.write('%14.4e ' %Kl)
            else:
                OF.write('%14.4f ' %Kl)
                
            if (k==5):
                k=0
                OF.write('\n') 
        OF.write('\n')
            
OF.close()
#============================ Plot
plt.figure(1)
plt.clf()
for g in range(0,G):
    plt.plot([Eg[g],Eg[g]],[plot_ymin,plot_ymax],color='k',linestyle='--')
plt.loglog(pyne_E,pyne_sigma_t,label=r'MT  1,total')
plt.loglog(pyne_E,pyne_sigma_s,label=r'MT  2,(z,elastic)')
plt.loglog(pyne_E,pyne_sigma_a,label=r'MT102,(z,$\gamma$)')
plt.loglog(StepE,StepXS_t,label='MT  1, MG',linestyle='--')
plt.loglog(StepE,StepXS_s,label='MT  2, MG',linestyle='--')
plt.loglog(StepE,StepXS_a,label='MT102, MG',linestyle='--')
plt.xlabel('Energy [MeV]')
plt.ylabel('Cross-section [b]')
plt.xlim(plot_xmin,plot_xmax)
plt.ylim(plot_ymin,plot_ymax)
plt.legend()
plt.show()