import os
import numpy as np 
import math
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy.optimize import fsolve
import module1 as m1 # module1 contains the function to convert r_star to r\n",
import cmath
from astropy.io import fits
import I_functions_class as Inp
from importlib import reload
import sys
import gc 
reload(Inp)

x = int(sys.argv[1])
omega_scale = float(sys.argv[2])


hdu_c = fits.open('/users/PCON0003/bowenchen12686/ondemand/data/sys/myjobs/electron_test/Constants.fits')
nu =  hdu_c[0].header['nu']
mu = hdu_c[0].header['mu']
lam = hdu_c[0].header['lam']
GC = hdu_c[0].header['GC']
c = hdu_c[0].header['c']
tol = hdu_c[0].header['tol']
direcPhoton = hdu_c[0].header['P_direc']
direcElectron = hdu_c[0].header['E_direc']
alpha = hdu_c[0].header['alpha']

M = hdu_c[0].header['M']
T = hdu_c[0].header['Temp']
omega=omega_scale*T
n=1

ks=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]
k_primes=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]
l_list = [1,2,3,4,5]
omega_index = round(omega*100*(8*np.pi*M))
print('omega_index',omega_index)

h = np.linspace(.01*T,20*T,2000)
h_index= round(h[x]*100*(8*np.pi*M))
h_prime1 = np.abs(h[x]-omega)
h_s = h[x]-omega



h_prime2 = h[x]+omega
h_prime1_index = round(np.abs(h[x]-omega)*100*(8*np.pi*M))
h_prime2_index = round((h[x]+omega)*100*(8*np.pi*M))

vals_A = np.zeros(8,dtype=complex)
vals_B = np.zeros(7,dtype=complex)
vals_C = np.zeros(7,dtype=complex)



# Initialize a list to temporarily store the data
#temp_data11[0] being empty
temp_data11 = [[] for _ in range(max(l_list) + 1)] #rpoints
temp_data12 = [[] for _ in range(max(l_list) + 1)] #F_points_up
temp_data13 = [[] for _ in range(max(l_list) + 1)] #z_points_up
temp_data14 = [[] for _ in range(max(l_list) + 1)] #F_points_in
temp_data15 = [[] for _ in range(max(l_list) + 1)] #z_points_in


for l in l_list:
    file_path = f"{direcPhoton}{l}ExtendedOmega.fits"
    with fits.open(file_path) as hduP:
        data11 = hduP[omega_index].data.field('rpoints_up')
        data12 = hduP[omega_index].data.field('F_points_up')
        data13 = hduP[omega_index].data.field('z_points_up')
        data14 = hduP[omega_index].data.field('F_points_in')
        data15 = hduP[omega_index].data.field('z_points_in')

        temp_data11[l].append(data11)
        temp_data12[l].append(data12)
        temp_data13[l].append(data13)
        temp_data14[l].append(data14)
        temp_data15[l].append(data15)

# Convert each list in temp_data to a NumPy array and store them in a NumPy array
r_points_gamma = np.array(temp_data11[1], dtype=float).squeeze()
psi_gammalomega_up = np.array([np.array(data, dtype=complex) for data in temp_data12 if len(data) > 0], dtype=object).squeeze(axis=1)
psi_gammalomega_prime_up = np.array([np.array(data, dtype=complex) for data in temp_data13 if len(data) > 0], dtype=object).squeeze(axis=1)
psi_gammalomega_in = np.array([np.array(data, dtype=complex) for data in temp_data14 if len(data) > 0], dtype=object).squeeze(axis=1)
psi_gammalomega_prime_in = np.array([np.array(data, dtype=complex) for data in temp_data15 if len(data) > 0], dtype=object).squeeze(axis=1)
rs = np.array([m1.r_star_to_r(x,M,tol) for x in r_points_gamma])
del temp_data11, temp_data12, temp_data13, temp_data14, temp_data15


temp_data21 = [[] for _ in range(11)] #F up+kh
temp_data22 = [[] for _ in range(11)] #F in+kh
temp_data23 = [[] for _ in range(11)] #G up+kh
temp_data24 = [[] for _ in range(11)] #G in+kh
temp_data25 = [[] for _ in range(11)] #R+k
temp_data26 = [[] for _ in range(11)] #T+k
temp_data27 = [[] for _ in range(11)] #F up+kh1
temp_data28 = [[] for _ in range(11)] #F in+kh1
temp_data29 = [[] for _ in range(11)] #G up+kh1
temp_data210 = [[] for _ in range(11)] #G in+kh1
temp_data211 = [[] for _ in range(11)] #F up+kh2
temp_data212 = [[] for _ in range(11)] #F in+kh2
temp_data213 = [[] for _ in range(11)] #G up+kh2
temp_data214 = [[] for _ in range(11)] #G in+kh2

for k in range(1, 11):
    file_path = f"{direcElectron}{k}ExtendedOmega.fits"
    with fits.open(file_path) as hdu:
        data21 = hdu[h_index].data.field('F_points_in')
        data22 = hdu[h_index].data.field('G_points_in')
        data23 = hdu[h_index].data.field('F_points_up')
        data24 = hdu[h_index].data.field('G_points_up')
        data25 = hdu[h_index].header['R']
        data26 = hdu[h_index].header['T']
        data27 = hdu[h_prime1_index].data.field('F_points_in')
        data28 = hdu[h_prime1_index].data.field('G_points_in')
        data29 = hdu[h_prime1_index].data.field('F_points_up')
        data210 = hdu[h_prime1_index].data.field('G_points_up')
        data211 = hdu[h_prime2_index].data.field('F_points_in')
        data212 = hdu[h_prime2_index].data.field('G_points_in')
        data213 = hdu[h_prime2_index].data.field('F_points_up')
        data214 = hdu[h_prime2_index].data.field('G_points_up')

        temp_data21[k].append(data21)
        temp_data22[k].append(data22)
        temp_data23[k].append(data23)
        temp_data24[k].append(data24)
        temp_data25[k].append(data25)
        temp_data26[k].append(data26)
        temp_data27[k].append(data27)
        temp_data28[k].append(data28)
        temp_data29[k].append(data29)
        temp_data210[k].append(data210)
        temp_data211[k].append(data211)
        temp_data212[k].append(data212)
        temp_data213[k].append(data213)
        temp_data214[k].append(data214)

F_points_xkh_in = np.array([np.array(data, dtype=complex) for data in temp_data21 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_in = np.array([np.array(data, dtype=complex) for data in temp_data22 if len(data) > 0], dtype=object).squeeze(axis=1)
F_points_xkh_up = np.array([np.array(data, dtype=complex) for data in temp_data23 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_up = np.array([np.array(data, dtype=complex) for data in temp_data24 if len(data) > 0], dtype=object).squeeze(axis=1)
Rk = np.array([np.array(data, dtype=complex) for data in temp_data25 if len(data) > 0], dtype=object).squeeze(axis=1)
Tk = np.array([np.array(data, dtype=complex) for data in temp_data26 if len(data) > 0], dtype=object).squeeze(axis=1)
F_points_xkh_prime1_in = np.array([np.array(data, dtype=complex) for data in temp_data27 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_prime1_in = np.array([np.array(data, dtype=complex) for data in temp_data28 if len(data) > 0], dtype=object).squeeze(axis=1)
F_points_xkh_prime1_up = np.array([np.array(data, dtype=complex) for data in temp_data29 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_prime1_up = np.array([np.array(data, dtype=complex) for data in temp_data210 if len(data) > 0], dtype=object).squeeze(axis=1)
F_points_xkh_prime2_in = np.array([np.array(data, dtype=complex) for data in temp_data211 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_prime2_in = np.array([np.array(data, dtype=complex) for data in temp_data212 if len(data) > 0], dtype=object).squeeze(axis=1)
F_points_xkh_prime2_up = np.array([np.array(data, dtype=complex) for data in temp_data213 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_prime2_up = np.array([np.array(data, dtype=complex) for data in temp_data214 if len(data) > 0], dtype=object).squeeze(axis=1)
del temp_data21, temp_data22, temp_data23, temp_data24,temp_data25,temp_data26,temp_data27,temp_data28,temp_data29,temp_data210,temp_data211,temp_data212,temp_data213,temp_data214


temp_data31 = [[] for _ in range(11)] #F up-kh
temp_data32 = [[] for _ in range(11)] #F in-kh
temp_data33 = [[] for _ in range(11)] #G up-khh
temp_data34 = [[] for _ in range(11)] #G in-khh
temp_data35 = [[] for _ in range(11)] #R-k
temp_data36 = [[] for _ in range(11)] #T-k
temp_data37 = [[] for _ in range(11)] #F up-kh1
temp_data38 = [[] for _ in range(11)] #F in-kh1
temp_data39 = [[] for _ in range(11)] #G up-kh1
temp_data310 = [[] for _ in range(11)] #G in-kh1
temp_data311 = [[] for _ in range(11)] #F up-kh2
temp_data312 = [[] for _ in range(11)] #F in-kh2
temp_data313 = [[] for _ in range(11)] #G up-kh2
temp_data314 = [[] for _ in range(11)] #G in-kh2


for k in range(1, 11):
    file_path = f"{direcElectron}min{k}ExtendedOmega.fits"
    with fits.open(file_path) as hdu:
        data31 = hdu[h_index].data.field('F_points_in')
        data32 = hdu[h_index].data.field('G_points_in')
        data33 = hdu[h_index].data.field('F_points_up')
        data34 = hdu[h_index].data.field('G_points_up')
        data35 = hdu[h_index].header['R']
        data36 = hdu[h_index].header['T']
        data37 = hdu[h_prime1_index].data.field('F_points_in')
        data38 = hdu[h_prime1_index].data.field('G_points_in')
        data39 = hdu[h_prime1_index].data.field('F_points_up')
        data310 = hdu[h_prime1_index].data.field('G_points_up')
        data311 = hdu[h_prime2_index].data.field('F_points_in')
        data312 = hdu[h_prime2_index].data.field('G_points_in')
        data313 = hdu[h_prime2_index].data.field('F_points_up')
        data314 = hdu[h_prime2_index].data.field('G_points_up')    

        temp_data31[k].append(data31)
        temp_data32[k].append(data32)
        temp_data33[k].append(data33)
        temp_data34[k].append(data34)
        temp_data35[k].append(data35)
        temp_data36[k].append(data36)
        temp_data37[k].append(data37)
        temp_data38[k].append(data38)
        temp_data39[k].append(data39)
        temp_data310[k].append(data310)
        temp_data311[k].append(data311)
        temp_data312[k].append(data312)
        temp_data313[k].append(data313)
        temp_data314[k].append(data314)
F_points_xminkh_in = np.array([np.array(data, dtype=complex) for data in temp_data31 if len(data) > 0], dtype=object)
G_points_xminkh_in = np.array([np.array(data, dtype=complex) for data in temp_data32 if len(data) > 0], dtype=object)
F_points_xminkh_up = np.array([np.array(data, dtype=complex) for data in temp_data33 if len(data) > 0], dtype=object)
G_points_xminkh_up = np.array([np.array(data, dtype=complex) for data in temp_data34 if len(data) > 0], dtype=object)
Rmink = np.array([np.array(data, dtype=complex) for data in temp_data35 if len(data) > 0], dtype=object)
Tmink = np.array([np.array(data, dtype=complex) for data in temp_data36 if len(data) > 0], dtype=object)
F_points_xminkh_prime1_in = np.array([np.array(data, dtype=complex) for data in temp_data37 if len(data) > 0], dtype=object)
G_points_xminkh_prime1_in = np.array([np.array(data, dtype=complex) for data in temp_data38 if len(data) > 0], dtype=object)
F_points_xminkh_prime1_up = np.array([np.array(data, dtype=complex) for data in temp_data39 if len(data) > 0], dtype=object)
G_points_xminkh_prime1_up = np.array([np.array(data, dtype=complex) for data in temp_data310 if len(data) > 0], dtype=object)
F_points_xminkh_prime2_in = np.array([np.array(data, dtype=complex) for data in temp_data311 if len(data) > 0], dtype=object)
G_points_xminkh_prime2_in = np.array([np.array(data, dtype=complex) for data in temp_data312 if len(data) > 0], dtype=object)
F_points_xminkh_prime2_up = np.array([np.array(data, dtype=complex) for data in temp_data313 if len(data) > 0], dtype=object)
G_points_xminkh_prime2_up = np.array([np.array(data, dtype=complex) for data in temp_data314 if len(data) > 0], dtype=object)
del temp_data31, temp_data32, temp_data33, temp_data34,temp_data35,temp_data36,temp_data37,temp_data38,temp_data39,temp_data310,temp_data311,temp_data312,temp_data313,temp_data314
###################################################

def FGsort(k,ksign,FGk,FGmink):
    if (k*ksign>0):
        return np.squeeze(FGk)
    else:
        return np.squeeze(FGmink)

def RTsort(k,RTk,RTmink):
    if (k>0):
        return np.array(RTk[k-1]).reshape(1)
    else:
        return np.array(RTmink[-k-1]).reshape(1)

def psort(k,k_prime,l):
    x = np.sign(k)*np.sign(k_prime)*(-1)**(k+k_prime+l)
    if ( x < 0):
        return 1
    else:
        return 0

def triangle_condition(j1, j2, j3):
    if (j1 + j2 >= j3) and (j1 + j3 >= j2) and (j2 + j3 >= j1):
        return 1
    else:
        return 0

def make_integrand():     
    integrand_A = 0
    integrand_B = 0
    integrand_C = 0

    array_zeros = np.zeros((20, 20, 5, 22))

    for k in range(len(ks)):
        j=(np.abs(ks[k]) - 1/2)
        
        for k_prime in range(len(k_primes)):
            j_prime=(np.abs(k_primes[k_prime]) - 1/2)
            for ll in range(len(l_list)):
                l = ll+1
                Tri = triangle_condition(j,j_prime,l)
                if Tri == 0:
                    continue
                pref = (2*j+1)/(2*l+1)
                if (h_s >0):
                    A1 = Inp.IfunctionsNoM(0,ks[k],1,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# I+-inupup h_prime1
                    A1_v = A1.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    A2 = Inp.IfunctionsNoM(1,ks[k],0,k_primes[k_prime],0,l,h[x],h_prime1,omega,M,n)# I+-upinin h_prime1
                    A2_v = A2.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    A3 = Inp.IfunctionsNoM(1,ks[k],1,k_primes[k_prime],0,l,h[x],h_prime1,omega,M,n)# I+-upupin h_prime1
                    A3_v = A3.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    A4 = Inp.IfunctionsNoM(1,ks[k],0,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# I+-upinup h_prime1
                    A4_v = A4.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]                
                    A5 = Inp.IfunctionsNoM(1,ks[k],1,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# A1**I+-upupup h_prime1
                    A5_v = A5.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]                
                    A6 = Inp.IfunctionsNoM(0,ks[k],0,k_primes[k_prime],0,l,h[x],h_prime1,omega,M,n)# A2*I+-ininin h_prime1*
                    A6_v = A6.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]                
                    A7 = Inp.IfunctionsNoM(0,ks[k],1,k_primes[k_prime],0,l,h[x],h_prime1,omega,M,n)# A3*I+-inupin h_prime1*
                    A7_v = A7.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]                
                    A8 = Inp.IfunctionsNoM(0,ks[k],0,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# A4*I+-ininup h_prime1*
                    A8_v = A8.IBarplusminsfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(k_primes[k_prime],1,F_points_xkh_prime1_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime1_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]                

                    array_zeros[k,k_prime,ll,0] = np.abs(RTsort(ks[k],Rk,Rmink))**2/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h_prime1)+1)*np.abs(A1_v)**2
                    array_zeros[k,k_prime,ll,1]  = -np.abs(RTsort(ks[k],Tk,Tmink))**2/(np.exp(8*np.pi*M*h_prime1)+1)*np.abs(A2_v)**2
                    array_zeros[k,k_prime,ll,2]  = -np.abs(RTsort(ks[k],Tk,Tmink))**2*(np.exp(8*np.pi*M*h_prime1))/(np.exp(8*np.pi*M*h_prime1)+1)/(np.exp(8*np.pi*M*h[x])+1)*np.abs(A3_v)**2
                    array_zeros[k,k_prime,ll,3]  =-np.abs(RTsort(ks[k],Tk,Tmink))**2*(np.exp(8*np.pi*M*omega))/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h[x])+1)*np.abs(A4_v)**2
                    array_zeros[k,k_prime,ll,4]  = np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h_prime1)+1)*np.conjugate(A1_v)*A5_v))[0].real
                    array_zeros[k,k_prime,ll,5]  = -np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))/(np.exp(8*np.pi*M*h_prime1)+1)*np.conjugate(A6_v)*A2_v))[0].real
                    array_zeros[k,k_prime,ll,6]  = -np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))*np.exp(8*np.pi*M*h_prime1)/(np.exp(8*np.pi*M*h_prime1)+1)/(np.exp(8*np.pi*M*h[x])+1)*np.conjugate(A7_v)*A3_v))[0].real
                    array_zeros[k,k_prime,ll,7]  = -np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))*np.exp(8*np.pi*M*omega)/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h[x])+1)*np.conjugate(A8_v)*A4_v))[0].real
                else: 
                    C1 = Inp.IfunctionsNoM(0,ks[k],0,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# I++ininup h_prime1
                    C1_v = C1.IBarplusplusfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(k_primes[k_prime],-1,F_points_xkh_prime1_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],-1,G_points_xkh_prime1_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    C2 = Inp.IfunctionsNoM(0,ks[k],1,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# I++inupup h_prime1
                    C2_v = C2.IBarplusplusfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(k_primes[k_prime],-1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],-1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    C3 = Inp.IfunctionsNoM(1,ks[k],1,k_primes[k_prime],0,l,h[x],h_prime1,omega,M,n)# I++upupin h_prime1
                    C3_v = C3.IBarplusplusfunc(ks[k],k_primes[k_prime],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(k_primes[k_prime],-1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],-1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    C4 = Inp.IfunctionsNoM(1,ks[k],0,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# I++upinup h_prime1
                    C4_v = C4.IBarplusplusfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(k_primes[k_prime],-1,F_points_xkh_prime1_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],-1,G_points_xkh_prime1_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_in[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    #C5 = C1* * C4 h_prime1
                    C6 = Inp.IfunctionsNoM(1,ks[k],1,k_primes[k_prime],1,l,h[x],h_prime1,omega,M,n)# C2**I++upupup h_prime1
                    C6_v = C6.IBarplusplusfunc(ks[k],k_primes[k_prime],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(k_primes[k_prime],-1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],-1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    C7 = Inp.IfunctionsNoM(0,ks[k],1,k_primes[k_prime],0,l,h[x],h_prime1,omega,M,n)# C3*I++inupin h_prime1
                    C7_v = C7.IBarplusplusfunc(ks[k],k_primes[k_prime],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h[x],h_prime1,omega,M,rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(k_primes[k_prime],-1,F_points_xkh_prime1_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],-1,G_points_xkh_prime1_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime1_up[abs(k_primes[k_prime])-1]))[psort(ks[k],k_primes[k_prime],l)]
                    
                    array_zeros[k,k_prime,ll,15]= np.abs(RTsort(ks[k],Rk,Rmink))**2/(np.exp(8*np.pi*M*omega)-1)*np.abs(C1_v)**2
                    array_zeros[k,k_prime,ll,16]= np.abs(RTsort(ks[k],Rk,Rmink))**2*np.exp(8*np.pi*M*h_prime1)/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h_prime1)+1)*np.abs(C2_v)**2
                    array_zeros[k,k_prime,ll,17]= -np.abs(RTsort(ks[k],Tk,Tmink))**2/(np.exp(8*np.pi*M*h_prime1)+1)/(np.exp(8*np.pi*M*h[x])+1)*np.abs(C3_v)**2
                    array_zeros[k,k_prime,ll,18]= np.abs(RTsort(ks[k],Tk,Tmink))**2*np.exp(8*np.pi*M*h[x])/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h[x])+1)*np.abs(C4_v)**2
                    array_zeros[k,k_prime,ll,19]= np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))*(2*np.exp(8*np.pi*M*h[x])+1)/(np.exp(8*np.pi*M*h[x])+1)/(np.exp(8*np.pi*M*omega)-1)*np.conjugate(C1_v)*C4_v))[0].real
                    array_zeros[k,k_prime,ll,20]= np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))*np.exp(8*np.pi*M*h_prime1)/(np.exp(8*np.pi*M*h_prime1)+1)/(np.exp(8*np.pi*M*omega)-1)*np.conjugate(C2_v)*C6_v))[0].real
                    array_zeros[k,k_prime,ll,21]=-np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))/(np.exp(8*np.pi*M*h_prime1)+1)/(np.exp(8*np.pi*M*h[x])+1)*np.conjugate(C2_v)*C6_v))[0].real


                B1 = Inp.IfunctionsNoM(1,k_primes[k_prime],0,ks[k],0,l,h_prime2,h[x],omega,M,n)# I+-upinin h_prime2
                B1_v = B1.IBarplusminsfunc(k_primes[k_prime],ks[k],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h_prime2,h[x],omega,M,rs,r_points_gamma,FGsort(k_primes[k_prime],1,F_points_xkh_prime2_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime2_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]))[psort(k_primes[k_prime],ks[k],l)]
                B2 = Inp.IfunctionsNoM(1,k_primes[k_prime],0,ks[k],1,l,h_prime2,h[x],omega,M,n)# I+-upinup h_prime2
                B2_v = B2.IBarplusminsfunc(k_primes[k_prime],ks[k],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h_prime2,h[x],omega,M,rs,r_points_gamma,FGsort(k_primes[k_prime],1,F_points_xkh_prime2_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime2_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]))[psort(k_primes[k_prime],ks[k],l)]                
                B3 = Inp.IfunctionsNoM(1,k_primes[k_prime],1,ks[k],0,l,h_prime2,h[x],omega,M,n)# I+-upupin h_prime2
                B3_v = B3.IBarplusminsfunc(k_primes[k_prime],ks[k],psi_gammalomega_in[ll],psi_gammalomega_prime_in[ll],l,h_prime2,h[x],omega,M,rs,r_points_gamma,FGsort(k_primes[k_prime],1,F_points_xkh_prime2_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime2_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]))[psort(k_primes[k_prime],ks[k],l)]                
                B4 = Inp.IfunctionsNoM(0,k_primes[k_prime],1,ks[k],1,l,h_prime2,h[x],omega,M,n)# I+-inupup h_prime2
                B4_v = B4.IBarplusminsfunc(k_primes[k_prime],ks[k],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h_prime2,h[x],omega,M,rs,r_points_gamma,FGsort(k_primes[k_prime],1,F_points_xkh_prime2_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime2_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime2_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime2_in[abs(k_primes[k_prime])-1]),FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]))[psort(k_primes[k_prime],ks[k],l)]                
                #B5 = B1*B3 h_prime2*
                B6 = Inp.IfunctionsNoM(0,k_primes[k_prime],0,ks[k],1,l,h_prime2,h[x],omega,M,n)# B4* * I+-ininup h_prime2
                B6_v = B6.IBarplusminsfunc(k_primes[k_prime],ks[k],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h_prime2,h[x],omega,M,rs,r_points_gamma,FGsort(k_primes[k_prime],1,F_points_xkh_prime2_in[abs(k_primes[k_prime])-1],F_points_xminkh_prime2_in[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime2_in[abs(k_primes[k_prime])-1],G_points_xminkh_prime2_in[abs(k_primes[k_prime])-1]),FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]))[psort(k_primes[k_prime],ks[k],l)]
                B7 = Inp.IfunctionsNoM(1,k_primes[k_prime],1,ks[k],1,l,h_prime2,h[x],omega,M,n)# B2*I+-upupup* h_prime2
                B7_v = B7.IBarplusminsfunc(k_primes[k_prime],ks[k],psi_gammalomega_up[ll],psi_gammalomega_prime_up[ll],l,h_prime2,h[x],omega,M,rs,r_points_gamma,FGsort(k_primes[k_prime],1,F_points_xkh_prime2_up[abs(k_primes[k_prime])-1],F_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(k_primes[k_prime],1,G_points_xkh_prime2_up[abs(k_primes[k_prime])-1],G_points_xminkh_prime2_up[abs(k_primes[k_prime])-1]),FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]))[psort(k_primes[k_prime],ks[k],l)]
                
                array_zeros[k,k_prime,ll,8] = np.abs(RTsort(ks[k],Rk,Rmink))**2/(np.exp(8*np.pi*M*h_prime2)+1)*np.abs(B1_v)**2
                array_zeros[k,k_prime,ll,9] = np.abs(RTsort(ks[k],Rk,Rmink))**2*np.exp(8*np.pi*M*omega)/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h_prime2)+1)*np.abs(B2_v)**2
                array_zeros[k,k_prime,ll,10] = np.abs(RTsort(ks[k],Tk,Tmink))**2*np.exp(8*np.pi*M*h[x])/(np.exp(8*np.pi*M*h_prime2)+1)/(np.exp(8*np.pi*M*h[x])+1)*np.abs(B3_v)**2
                array_zeros[k,k_prime,ll,11] = -np.abs(RTsort(ks[k],Tk,Tmink))**2/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h[x])+1)*np.abs(B4_v)**2
                array_zeros[k,k_prime,ll,12] = np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))*(2*np.exp(8*np.pi*M*h[x])+1)/(np.exp(8*np.pi*M*h_prime2)+1)/(np.exp(8*np.pi*M*h[x])+1)*np.conjugate(B3_v)*B1_v))[0].real
                array_zeros[k,k_prime,ll,13] = -np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h[x])+1)*np.conjugate(B4_v)*B6_v))[0].real
                array_zeros[k,k_prime,ll,14] = np.array(((np.conjugate(RTsort(ks[k],Rk,Rmink))*RTsort(ks[k],Tk,Tmink))*np.exp(8*np.pi*M*omega)/(np.exp(8*np.pi*M*omega)-1)/(np.exp(8*np.pi*M*h_prime2)+1)*np.conjugate(B2_v)*B7_v))[0].real



                integrand_A += pref*np.sum(array_zeros[k, k_prime, ll, 0:8])
                integrand_B += pref*np.sum(array_zeros[k, k_prime, ll, 8:15])
                integrand_C += pref*np.sum(array_zeros[k, k_prime, ll, 15:])

        
    return (integrand_A*alpha/(4*np.pi*np.pi),integrand_B*alpha/(4*np.pi*np.pi),integrand_C*alpha/(4*np.pi*np.pi),array_zeros*alpha/(4*np.pi*np.pi))

integrand_A,integrand_B,integrand_C,array_zeros = make_integrand()

# Save the array and metadata for this energy scale
hdu = fits.PrimaryHDU(data=array_zeros)

# Add metadata to the header
hdu.header['INTEG_A'] = integrand_A
hdu.header['INTEG_B'] = integrand_B
hdu.header['INTEG_C'] = integrand_C
hdu.header['x'] = x
hdu.header['OMEGA'] =omega_scale 

# Save to a FITS file (unique filename for each energy scale)
hdu.writeto(f'output_x{x}_omega{omega_scale}.fits', overwrite=True)

del r_points_gamma,psi_gammalomega_up,psi_gammalomega_prime_up,psi_gammalomega_in,psi_gammalomega_prime_in,rs
del F_points_xkh_in,G_points_xkh_in,F_points_xkh_up,G_points_xkh_up,Rk,Tk,F_points_xkh_prime1_in,G_points_xkh_prime1_in,F_points_xkh_prime1_up,G_points_xkh_prime1_up,F_points_xkh_prime2_in,G_points_xkh_prime2_in,F_points_xkh_prime2_up,G_points_xkh_prime2_up
del F_points_xminkh_in,G_points_xminkh_in,F_points_xminkh_up,G_points_xminkh_up,Rmink,Tmink,F_points_xminkh_prime1_in,G_points_xminkh_prime1_in,F_points_xminkh_prime1_up,G_points_xminkh_prime1_up,F_points_xminkh_prime2_in,G_points_xminkh_prime2_in,F_points_xminkh_prime2_up,G_points_xminkh_prime2_up

gc.collect()