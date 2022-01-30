# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:56:01 2022

@author: rapha
"""

from math import tan,atan2,cos,sin
import numpy as np
import sys 
import matplotlib.pyplot as plt
from scipy import fftpack
from copy import *
from pylab import *

#%%

PI=3.1415926535898;


def computeCk(U_inf,chord):     
        Omega = 450
        k = Omega*chord/U_inf
        Z0 = complex(1 , (0.0455/k))
        Z1 = complex(1 , (0.3/k))
        Z3 = 1 - (0.165/Z0) - (0.335/Z1) 
        Ck = ((Z3.real)**2 + (Z3.imag)**2)**0.5
        return Ck


#--------------------------------------------------------------------------------

def ClaCmaAilerons(d,a,b,s,m,n):
    """d diamètre à l'ogive, a emplenture ailerons, b saumon, s envergure, 
    m flèche(des ailerons),n nombre ailerons (3 ou 4 en config finale)"""

## la longueur de référence est l'emplanture et la surface celle de la section du tube

    #ailerons
    e = m+b-a
    z=((m-a/2+b/2)**2+s**2)**1/2
    ail=[e+(m*(a+2*b)/(3*(a+b))+1/6*(a+b-(a*b/(a+b)))),(4*n*(s/d)**2)/(1+(1+((2*z)/(a+b))**2)**1/2)]
    
    if n == 4 :
        ail_interference=[ail[0],ail[1]*(1+d/(2*s+d))]
    else :
        ail_interference = ail
    
    cdm = (m/2*(2*m+a)+(m+b-a/2)*2*a+b)/3*(m+a)

    Cla = ail_interference[1]
    Cma = ail_interference[1]*(cdm-ail_interference[0])/a

    return [Cla, Cma]

""" """
def computeLift(U_inf,Surf,chord,VecH,VecTeta,dCLa):

        """
        Fonction mise à jour avec dCLa en suivant Amandolese et al.
        """        
        #rho = np.sqrt(2 * q / (U_inf * U_inf)) --> constante utilisé avant pour effectuer le dCL/dAlpha
        b = chord / 2
        a = S0 / (Mass * b)
        #OLD --> FA = rho * PI * b * b * (VecTeta[1] * U_inf + VecH[2] - a * b * VecTeta[2])
        Cz = PI * b * b * (VecTeta[1]  / U_inf + VecH[2] / U_inf**2 - a * b * VecTeta[2] / U_inf**2)

        #OLD --> FA += 2 * PI * computeCk(U_inf, chord) * rho * U_inf * b * (U_inf * VecTeta[0] + VecH[1] + (0.5 - a) * b * VecTeta[1])
        Cz +=  dCLa* computeCk(U_inf, chord) * b * (VecTeta[0] + VecH[1] / U_inf + (0.5 - a) * b * VecTeta[1] / U_inf)

        #------------------
        
        return Cz
    

    
def computeMoment(U_inf,Surf,chord,VecH,VecTeta,dCMa):

        # Pure steady model
        #------------------
        #CM = 2*PI*(VecTeta[0] + VecH[1]/U_inf)/4.0; #to be modified for part 2 and part 3
        #rho = np.sqrt(2 * q / (U_inf * U_inf))
        
        b = chord / 2
        a = S0 / (Mass * b)
        #OLD --> MA = rho * PI * b * b * (a * b * VecH[2] + U_inf * b * (a - 0.5) * VecTeta[1] - b * b * (0.125 + a * a) * VecTeta[2])
        Cm = (PI/2) * b**2 * (a * b * VecH[2] / U_inf**2 - b**2 * (0.125 + a**2) * VecTeta[2] / U_inf**2 + b * (a - 0.5) * VecTeta[1] / U_inf)
        #OLD --> MA += 2 * PI * computeCk(U_inf, chord) * rho * U_inf * b * b * (U_inf * VecTeta[0] + VecH[1] + (0.5 - a) * b * VecTeta[1])
        Cm += dCMa*computeCk(U_inf, chord) *b**2 * ( VecTeta[0] + VecH[1] / U_inf + (0.5 - a) * b * VecTeta[1] / U_inf)
        #------------------

        return Cm

def UpdateVP(Mass,CouplingTerm,DissTerm,StiffTerm,AeroTerm,XI):
   A = [];
   A.append(0.0);
   A.append(0.0);
   A[0] = XI[1];
   A[1] = (1/Mass)*(AeroTerm - CouplingTerm - StiffTerm - DissTerm)
   return(A)

def IntegrateRK4(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, Vec1, DT):

        X0 = [];
        X0.append(Vec1[0]);
        X0.append(Vec1[1]);
        X1 = [];
        X1.append(0.0);
        X1.append(0.0);
        X2 = [];
        X2.append(0.0);
        X2.append(0.0);
        X3 = [];
        X3.append(0.0);
        X3.append(0.0);
        A0 = [];
        A0.append(0.0);
        A0.append(0.0);
        A1 = [];
        A1.append(0.0);
        A1.append(0.0);
        A2 = [];
        A2.append(0.0);
        A2.append(0.0);
        A3 = [];
        A3.append(0.0);
        A3.append(0.0);

        A0[0] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X0)[0];
        A0[1] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X0)[1];
        X1[0] = X0[0] + A0[0]*DT/2.0;
        X1[1] = X0[1] + A0[1]*DT/2.0;

        A1[0] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X1)[0];
        A1[1] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X1)[1];
        X2[0] = X0[0] + A1[0]*DT/2.0;
        X2[1] = X0[1] + A2[1]*DT/2.0;

        A2[0] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X2)[0];
        A2[1] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X2)[1];
        X3[0] = X0[0] + A2[0]*DT/1.0;
        X3[1] = X0[1] + A2[1]*DT/1.0;

        A3[0] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X3)[0];
        A3[1] = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X3)[1];

        Vec1[0] = Vec1[0] + (DT/6.0)*(A0[0]+2.0*A1[0]+2.0*A2[0]+A3[0]);
        Vec1[1] = Vec1[1] + (DT/6.0)*(A0[1]+2.0*A1[1]+2.0*A2[1]+A3[1]);
        Vec1[2] = (1./6.)*(A0[1]+2.0*A1[1]+2.0*A2[1]+A3[1]) 

        return(Vec1)

#%%

IsH = 1 #Activate time marching for plunging at IsH=1
IsTeta = 1 #Same with IsTheta
U_inf = 50
U_kmh = U_inf * 60 * 60 / 1000

# Fluid model
#------------
#rho = 1.205 OLD
chord = 0.22
span = 0.3
Surf = chord*span/2
#q = (1./2.)*rho*(U_inf**2)
Teta0 = 0.0001
h0 = 0.001
xCG = 0.08 * chord #8% basé sur BE1_flutter


# Parametres géometriques

""" d diamètre à l'ogive
    a emplenture ailerons
    b saumon
    s envergure, 
    m flèche(des ailerons)
    n nombre ailerons (3 ou 4 en config finale),
    --> En metres
"""
## Données SERA-IV Nov18

d = 0.156
a = 0.35
b = 0.09
s = 0.3
m = 0.397
n = 3


# Structure model
#----------------
omega_h = 570
omega_a = 412
         
StiffH = 35500    #Redefinition sur catia 35500
StiffTeta = 1216   #Redefinition sur catia 1216
DissH = 0.0 #Redefinition sur catia
DissTeta = 0.0    #Redefinition sur catia
#Diss_ratio = DissH/DissTeta

Mass = StiffH / omega_h**2



# Others derived quantities
#--------------------------
S0 = Mass*xCG
Inertia = StiffH / omega_a**2
PI = 3.1415926535898

# Numerical parameters
#---------------------
# Number of time steps 
NbIt = 400000
MaxTime = 5.
DT = MaxTime/NbIt





#%%

#---------------------------------------------------------
# Here we initialize both vectors of coordinates
#---------------------------------------------------------
VecH = []
VecTeta = []
VecH.append(h0)
VecH.append(0.0)
VecH.append(0.0)
VecTeta.append(Teta0)
VecTeta.append(0.0)
VecTeta.append(0.0)

#---------------------------------------------------------
# We create variables to monitore both
# displacement and velocity for each of
# the two coordinates
#---------------------------------------------------------
DispH = []
DispTeta = []
VelH = []
VelTeta = []
VecTime = []

#---------------------------------------------------------
#Calculate dCMa and dCLa
dCLa, dCMa = ClaCmaAilerons(d,a,b,s,m,n)

Lifting = []


study_list = []



data_dict = {
        'Dh'                : 0.0,
        'Da'                : 0.0,
        'Dh/Da'             : 0.0,
        'FreqH'             : 0.0,
        'FreqTeta'          : 0.0,
        'FreqH/FreqTeta'    : 0.0
    }


# Loop on time steps
for CurrTime in range(0,NbIt):
    
    Time = CurrTime*DT
    VecHNEW = []
    VecTetaNEW = []
    #----------------------------------------------
    # Update of the profile positions (h, teta) 
    # displacement velocity (dot(h), dot(teta)) 
    # and acceleration (ddot(h), ddot(teta)) 
    #----------------------------------------------
    
    Lift = 1.0*computeLift(U_inf,Surf,chord,VecH,VecTeta,dCLa)
    CouplingTerm = IsTeta*S0*VecTeta[2]
    StiffTerm = StiffH*VecH[0]
    DissTerm = DissH*VecH[1]
    if (IsH):
     VecHNEW = IntegrateRK4(Mass, CouplingTerm, DissTerm, StiffTerm, Lift, VecH, DT)
     
     
     VecH = VecHNEW
     #print("working")
    else:
     VecHNEW = VecH
     VecH = VecHNEW
     
    Lifting.append(Lift)
     
    Moment = 1.0*computeMoment(U_inf,Surf,chord,VecH,VecTeta,dCMa)
    CouplingTerm = IsH*S0*VecH[2]
    StiffTerm = StiffTeta*(VecTeta[0])
    DissTerm = DissTeta*VecTeta[1]
    if (IsTeta):
      VecTetaNEW = IntegrateRK4(Inertia, CouplingTerm, DissTerm, StiffTerm, Moment, VecTeta, DT)
      VecTeta = VecTetaNEW
    else:
      VecTetaNEW = VecTeta
      VecTeta = VecTetaNEW

    # Data to save for post-processing
    #---------------------------------
    DispH.append(VecH[0])
    VelH.append(VecH[1])
    DispTeta.append(VecTeta[0])
    VelTeta.append(VecTeta[1])
    VecTime.append(Time)


# Estimation of frequencies
#--------------------------
DispH2 = []
DispTeta2 = []
Ind1 = int(NbIt/2)-1
DispH2 = DispH[Ind1:-1]
DispTeta2 = DispTeta[Ind1:-1]
#
FFT1 = abs(fftpack.fft(DispH2))
FFTH = FFT1.tolist()
FFT2 = abs(fftpack.fft(DispTeta2))
FFTTeta = FFT2.tolist()
#
FFTH[0]=0.
FFTTeta[0]=0.
FreqH = 0.0
FreqTeta = 0.0
FreqH = (FFTH.index(max(FFTH[0:(np.size(FFTH)-1)//2])))/(MaxTime/2.)
FreqTeta = (FFTTeta.index(max(FFTTeta[0:(np.size(FFTTeta)-1)//2])))/(MaxTime/2.)
print("U_inf   = ",U_inf, "m/s  ->", U_kmh, "km/h")
print("Dh/Da  =", round(DissRatio, 4))
print("FreqH   = ",FreqH)
print("FreqTeta= ",FreqTeta)
print("FreqH/FreqTeta=  ", round(FreqH/FreqTeta, 4))
#
# Graphic functions
#------------------
plt.figure(figsize=(10,7))

plt.title('Aileron Flutter Test')
plt.plot(VecTime,DispH,"b-", label="h")
plt.xlabel('Time t')
plt.ylabel('Displacement h')
plt.legend(loc='upper left')
plt2 = plt.twinx()
plt2.plot(VecTime,DispTeta,"r-", label="Teta")
plt2.set_ylabel('Displacement Teta')
plt2.legend()
plt.show()

data_dict['Dh'] = DissH
data_dict['Da'] = DissTeta
data_dict['Dh/Da'] = DissRatio
data_dict['FreqH'] = FreqH
data_dict['FreqTeta'] = FreqTeta
data_dict['FreqH/FreqTeta'] = round(FreqH/FreqTeta, 4)

study_list.append(data_dict)

DissRatio += 2

