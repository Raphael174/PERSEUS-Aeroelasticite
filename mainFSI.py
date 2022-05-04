# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:44:26 2022

@author: Perseus - Raphaël AUBRY, Eliot MARTIN, Adrien SIMONNET

########################## MODELE ADIMENSIONE ################################

Parametres d'entrées
--------------------

Géometrie de l'aile     

c = corde
b = semi_corde
a =  semi-corde au centre elastique (%)
r_alpha = rayon de gyration adimensioné
x_CG = distance de CG à CE m
x_alpha = distance adimensioné de CG à CE 

Dans X. Amandolese et al. (2013), x_CG = a


Structure               
omega_teta = fréquence propre torsion rad/s
omega_h = fréquence propre flexion    rad/s
zeta_alpha = amortissement structure torsion
zeta_h = amortissement structure flexion


Aero 
rho = masse volumique de l'air kg/m^3
M = Mach
q_inf = pression dynamique (ici, rho*b*U_inf**2*Surf/2)
t_0 = ratio b/U_inf (provenent du rapport de Claudia)
k =     k = omega_teta*b/U_inf constante de theodorsen
dCl = tangente Cl/alpha, 2*PI pour une plaque 2D
dCm = tangente Cm/alpha

Valeurs intiales en h et teta
-avec a = une petite valeur (~ teta = 1° ou h = 1cm)
-On remarque que lorsque l'on a h0 = 0 et teta0 = a, une oscillation forte de la
 part de teta est observé au début, et des frequences très différentes (omega_teta ~ 10*omega_h)
-On remarque aussi que lorsque h0 = a et teta0 = 0, les frequences d'oscillations restent égales
 et les deux mouvement se suivent très clairement
== Il est donc préférable d'utilsé une valeur non-zero pour teta0 et zero pour h0
    avec teta0 = x*PI/180 (x = 1, 2, 3 °)

"""

import matplotlib.pyplot as plt
import pandas as pd
from timeit import default_timer as timer
import numpy as np


from data_process_functions import computePlotDamping # takes hDisp : list, tetaDisp: list , timeList: list
from data_process_functions import computeFrequencies # takes Nb, posH, posTeta, MaxTime
from integration_functions import IntegrateRK4       # takes Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, Vec1, DT
from aero_functions import computeLift # takes omega_teta, rho, a, b, c, U_inf, Surf, half_c, VecH, VecTeta
from aero_functions import computeMoment #takes omega_teta, rho, a, b, c, U_inf, Surf, half_c, VecH, VecTeta

PI = 3.1415926535898

plot_if = True #Change to True to display data
to_excel = False #send frequency and velocity data to excel

plot_damping = False

plot_frequency = True

print("Simulation start ------------ \n")
start = timer()

#%% DONNEES D'ENTREE

comp = 1; #mettre à 1 pour prendre en compte le compressible
IsH = 1  # Activate time marching for plunging at IsH=1
IsTeta = 1 # Same with IsTheta


# Structure model
# ---------------------
m = 0.304
StiffH = 595.6
StiffTeta = 0.149

omega_h = 44.26
zeta_h = .0 #5.38e-2
DissH = .0
omega_teta = 56.55
zeta_teta = .0 #7.91e-5
DissTeta = .0

# Geometry
# ---------------------

c = 0.035
a = 2.8e-3  #
b = c/2
x_alpha = -0.159 # For coupling
S0 = -1.0*8.5e-4
r_alpha = 0.707 #radius of gyration
s = 0.225
Surf = c*s

# Aero
# ---------------------
rho = 1.205
T_air = (15 + 273.15)
v_sound = np.sqrt(1.4 * 8.314e3 * T_air)
dCl =  2*PI #6.2
dCm = ((b-2.8e-3)/c-0.25)*2*PI #0.91


# Integration parameters
# ---------------------
DT = 5 / 400000
MaxTime = 2.5
NbIt = int(MaxTime/DT) #Put here for MaxTime study
# NbIt = 400000
# MaxTime = 5.
# DT = MaxTime/NbIt

# Initial positions
# ---------------------
x = 1
z = 0
Teta0 = 0.017453 #x*PI/180 # PI/180 = 1°
h0 = 0.001 #z*.01


#%% MAIN

Fh = []
Fa = []
V = []

frequency_h_ratio = []
frequency_teta_ratio = []
damping_h = []
damping_teta = []

velocity = []

U_inf = 0.5
run = 1
#time_array = []

while U_inf <=15:   
                  
    # initialize vectors 
    VecH = [h0, 0.0, 0.0]
    VecTeta = [Teta0, 0.0, 0.0]
    
    # data storage
    DispH = []
    DispTeta = []
    VelH = []
    VelTeta = []
    VecTime = []
        
    # Loop on time steps
    for CurrTime in range(0, NbIt):
        """ Code that is commented out is Claudia's version of the terms"""
    
        Time = CurrTime * DT
        VecHNEW = []
        VecTetaNEW = []
        
        ################# H
        
        Lift =  computeLift(omega_teta, rho, dCl, a, b, c, U_inf, Surf, VecH, VecTeta, comp) #1/(m*b) * computeLift(omega_teta, rho, dCl, a, b, c, U_inf, Surf, VecH, VecTeta)
        CouplingTerm = IsTeta*S0*VecTeta[2] #IsTeta*x_alpha* VecTeta[2] 
        StiffTerm =  StiffH*VecH[0] #omega_h**2 * VecH[0]
        DissTerm =  DissH*VecH[1] #2*zeta_h*omega_h**2 * VecH[1]
        MassTerm =  m
        
        if IsH:
            VecHNEW = IntegrateRK4(MassTerm, CouplingTerm, DissTerm, StiffTerm, Lift, VecH, DT)
            VecH = VecHNEW
        else:
            VecHNEW = VecH
            VecH = VecHNEW
    
        ################# TETA
        
        Moment = 1.0*computeMoment(omega_teta, rho, dCm, a, b, c, U_inf, Surf, VecH, VecTeta, comp) # 1/(m*b**2) * computeMoment(omega_teta, rho, dCm, a, b, c, U_inf, Surf, VecH, VecTeta)
        CouplingTerm =  IsH*S0*VecH[2] #sH * x_alpha * 1/(m*b) * VecH[2]
        StiffTerm = StiffTeta*(VecTeta[0]) # omega_teta**2 * r_alpha**2 * VecTeta[0]
        DissTerm =  DissTeta*VecTeta[1] #2*zeta_teta * r_alpha**2 * omega_teta**2 * VecTeta[1]
        InertiaTerm =  4.66e-5 #1/r_alpha**2
        
        if IsTeta:
            VecTetaNEW = IntegrateRK4(InertiaTerm, CouplingTerm, DissTerm, StiffTerm, Moment, VecTeta, DT)
            VecTeta = VecTetaNEW
        else:
            VecTetaNEW = VecTeta
            VecTeta = VecTetaNEW
    
        # data storage
        DispH.append(VecH[0]) #round(VecH[0]*1000, 1)
        VelH.append(VecH[1])
        DispTeta.append(VecTeta[0]) #round(VecTeta[0]*180/PI, 1)
        VelTeta.append(VecTeta[1])
        VecTime.append(Time)
        
        # Barre de chargement
        if CurrTime % (NbIt // 10) == 0:
            print(f"{CurrTime // (NbIt // 100)} %")
    
        
    #########################################################################
    ############################# DATA PROCESSING ###########################
    # calcule de l'amortissement en h et teta à chaque simulation
 
    dampingOutput_Data                  = computePlotDamping(DispH, DispTeta, VecTime)
    #avg_damping_h, var_damping_h        = .0
    avg_damping_h, var_damping_h        = dampingOutput_Data[0], dampingOutput_Data[1]
    #avg_damping_teta, var_damping_teta  = .0
    avg_damping_teta, var_damping_teta  = dampingOutput_Data[2], dampingOutput_Data[3]
    
    # verification de la variance des amortissement
    if var_damping_h > 0.05:
        print('h damping variance is above 0.05')
    if var_damping_teta > 0.05:
        print('teta damping variance is above 0.05')
    

    # Estimation of frequencies
    FreqH, FreqTeta = computeFrequencies(NbIt, DispH, DispTeta, MaxTime)

    # Storing current parameters
    current_frequency_h = round(FreqH, 2)
    current_frequency_teta = round(FreqTeta, 2)
    
    # adding parameters to lists
    frequency_h_ratio.append(current_frequency_h)
    frequency_teta_ratio.append(current_frequency_teta)
    damping_h.append(avg_damping_h)
    damping_teta.append(avg_damping_teta)
    velocity.append(U_inf)
    
    #time_array.append(MaxTime) #Append Max Time
    
    end = timer() #fin de chrono
    
    # debugging
    print("------------------------")
    print('Run °', run) 
    print("Running time in s= ",round(end-start, 2), "\n")
    print("MaxTime selected in s = ", MaxTime)
    print("\n")
    print('FLEXION')
    print("Frequency in h = ", current_frequency_h)
    print("Damping in h = ", avg_damping_h)
    print("\n")
    print('TORSION')
    print("Frequency in teta = ", current_frequency_teta)
    print("Damping in teta = ", avg_damping_teta)
    print("\n")
    print("U_inf = ", U_inf)
    print("------------------------ \n\n")
    

    
    #####################################################################
    # Graphic functions
    # ------------------

        
    if plot_if: 
        plt.figure(figsize=(8, 6))
        plt.title('@U_inf = %1.3f' %U_inf)
        plt.plot(VecTime, DispH, color = 'blue', label=("(gloabl) omega_h = %1.1f Hz" %FreqH))
        plt.xlabel('Time (s)')
        plt.ylabel('Heave (m)')
        plt.legend(loc='upper left', frameon=True)
        plt2 = plt.twinx()
        plt2.plot(VecTime, DispTeta, color = 'green', label=("(global) omega_teta: %1.1f Hz" %FreqTeta))
        plt2.set_ylabel('Pitch (°)')
        plt2.legend(loc='upper right', frameon=True)
        plt.show()
    
    Fh.append(FreqH)
    Fa.append(FreqTeta)
    V.append(U_inf)
    
    run += 1
    U_inf += 0.5 ########## LOOP OVER
   
#%%
     
if plot_damping:
    plt.figure(figsize=(8, 6))
    plt.title('Growth rate vs. Air Speed')
    plt.plot(velocity, damping_h, "b-", label=("h"))
    plt.xlabel('Air Speed (m/s)')
    plt.ylabel('Damping_h')
    plt.legend(loc='upper left', frameon=True)
    plt2 = plt.twinx()
    plt2.plot(velocity, damping_teta, "r-", label=("teta"))
    plt2.set_ylabel('damping_teta')
    plt2.legend(loc='upper right', frameon=True)
    plt.show()
    
if plot_frequency:
    plt.figure(figsize=(8, 6))
    plt.title('Frequency vs. Air Speed')
    plt.plot(V, Fh, color = 'blue', label=("Heave frequency"))
    plt.xlabel('U_inf m/s')
    plt.ylabel('omega_h (Hz)')
    plt.legend(loc='upper left', frameon=True)
    plt2 = plt.twinx()
    plt2.plot(V, Fa, color = 'green', label=("Pitch frequency"))
    plt2.set_ylabel('omega_teta (Hz)')
    plt2.legend(loc='upper right', frameon=True)
    plt.show()

    
#%%

if to_excel:
    col1 = "FreqH"
    col2 = "FreqTeta"
    col3 = "Velocity"
    data = pd.DataFrame({col1:frequency_h_ratio, col2:frequency_teta_ratio, col3:velocity})
    data.to_excel('Results_2.5s_Uinf_1-15_FullFFT_Amendolese_usingClaudia.xlsx', sheet_name='sheet1', index=False)
    
