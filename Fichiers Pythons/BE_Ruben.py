# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 14:49:22 2022

@author: rapha
"""

from math import tan,atan2,cos,sin
import numpy as np
import sys 
import matplotlib.pyplot as plt
from scipy import fftpack
from copy import *
from pylab import *

# definition of problem parameters

PI=3.1415926535898;


#For the part 3 : the function $computeCk$ is a first assessment of the value of C_k=R(t)/(2Q(t)).
#-------------------------------------------------------------------------------------------------------------------------------------

def computeCk(U_inf,chord):     
        Omega = 56
        k = Omega*chord/U_inf
        Z0 = complex(1 , (0.0455/k))
        Z1 = complex(1 , (0.3/k))
        Z3 = 1 - (0.165/Z0) - (0.335/Z1) 
        Ck = ((Z3.real)**2 + (Z3.imag)**2)**0.5
        return Ck

#The following function is designed to compute the lift F_A 
#--------------------------------------------------------------------------------

def computeLift(q,U_inf,Surf,chord,VecH,VecTeta):

        # Pure steady model
        #------------------
        #CZ = -2.0*PI*(VecTeta[0] + VecH[1]/U_inf) #to be modified for part 2 and part 3
        rho = np.sqrt(2 * q / (U_inf * U_inf))
        b = chord / 2
        a = S0 / (Mass * b)
        FA = rho * PI * b * b * (VecTeta[1] * U_inf + VecH[2] - a * b * VecTeta[2])
        FA += 2 * PI * computeCk(U_inf, chord) * rho * U_inf * b * (U_inf * VecTeta[0] + VecH[1] + (0.5 - a) * b * VecTeta[1])
        #------------------
        
        return(FA)
    

    
def computeMoment(q,U_inf,Surf,chord,VecH,VecTeta):

        # Pure steady model
        #------------------
        #CM = 2*PI*(VecTeta[0] + VecH[1]/U_inf)/4.0; #to be modified for part 2 and part 3
        rho = np.sqrt(2 * q / (U_inf * U_inf))
        b = chord / 2
        a = S0 / (Mass * b)
        MA = rho * PI * b * b * (a * b * VecH[2] + U_inf * b * (a - 0.5) * VecTeta[1] - b * b * (0.125 + a * a) * VecTeta[2])
        MA += 2 * PI * computeCk(U_inf, chord) * rho * U_inf * b * b * (U_inf * VecTeta[0] + VecH[1] + (0.5 - a) * b * VecTeta[1])
        #------------------

        return(MA)
    
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
    
IsH = 1 #Activate time marching for plunging at IsH=1
IsTeta = 1 #Same with IsTheta

# Fluid model
#------------
rho = 1.205
chord = 0.035
span = 0.23
Surf = chord*span/2
#q = (1./2.)*rho*(U_inf**2)
Teta0 = 0.017453
h0 = 0.001

# Structure model
#----------------
Mass = 0.304
StiffH = 595.6
StiffTeta = 0.149
DissH = 0.0538
DissTeta = 0.0000791

# Others derived quantities
#--------------------------
S0 = 0.00085
Inertia = 4.66e-5
PI = 3.1415926535898

# Numerical parameters
#---------------------
# Number of time steps 
NbIt = 400000
MaxTime = 5.
DT = MaxTime/NbIt

U_inf = 10

Fh = []
Fa = []
V = []

run = 1

while U_inf <= 100:   
    
    q = (1./2.)*rho*(U_inf**2)
    
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
        Lift = 1.0*computeLift(q,U_inf,Surf,chord,VecH,VecTeta)
        CouplingTerm = IsTeta*S0*VecTeta[2]
        StiffTerm = StiffH*VecH[0]
        DissTerm = DissH*VecH[1]
        if (IsH):
         VecHNEW = IntegrateRK4(Mass, CouplingTerm, DissTerm, StiffTerm, Lift, VecH, DT)
         VecH = VecHNEW
        else:
         VecHNEW = VecH
         VecH = VecHNEW
         
        Moment = 1.0*computeMoment(q,U_inf,Surf,chord,VecH,VecTeta)
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
    print("U_inf   = ",U_inf)
    print("FreqH   = ",FreqH)
    print("FreqTeta= ",FreqTeta)
    #
    # Graphic functions
    #------------------
    plt.figure(figsize=(12,8))
    plt.title('FSI-Airfoil')
    plt.plot(VecTime,DispH,"b-", label="Coord. h")
    plt.xlabel('Time t')
    plt.ylabel('Displacement h')
    plt.legend()
    plt2 = plt.twinx()
    plt2.plot(VecTime,DispTeta,"r-", label="Coord. Teta")
    plt2.set_ylabel('Displacement Teta')
    plt2.legend()
    plt.show()
    
    Fh.append(FreqH)
    Fa.append(FreqTeta)
    V.append(U_inf)
    
    U_inf += 10
    
    print('Run', run)
    
    run += 1
    
    