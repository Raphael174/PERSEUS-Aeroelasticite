# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 18:35:39 2022

@author: Perseus - Raphaël AUBRY, Eliot MARTIN, Adrien SIMONNET
"""

import numpy as np
from scipy import fftpack

PI = 3.1415926535898


def computePlotDamping(hDisp : list, tetaDisp: list , timeList: list):
    """ Dernier ajout au BE, le 20/04/2022 par Raphaël
        Le but est de calculer l'amortissement observé sur les courbes 
        oscillatoires de flexion et torsion pour les comparer à des études """
    peaks_h     = []    #calcul des point max en h
    peaks_teta  = []    #calcul des point max en h  
    n = 5               # nombre de point considéré pour trouver une valeur max
    b = 10
        
    for i in range(b*n, len(hDisp)-b*n):
        if max(hDisp[i-n:i+n]) == hDisp[i]: #le tri des valeurs max
            peaks_h.append(i)
        # if i >=10: #pour eviter de passer en revue trop de données
        #     break
    for i in range(n, len(tetaDisp)-n):
        if max(tetaDisp[i-n:i+n]) == tetaDisp[i]:
            peaks_teta.append(i)
        # if i >=10:
        #     break
            
    hMax = []
    tetaMax = []
    TimeMaxH = []
    TimeMaxTeta = [] 
    """Creating two separate loops since h and teta may not 
    necessarily have the same number of peaks in a given simulation
    Also, multiple if-elseif statement to only take +ive consecutive values """
    for i in range(1, len(peaks_h)):
        if len(hMax) >= 5: 
            break #first if-statement to decide if there are enough data points
        elif hDisp[peaks_h[i]] > 0:
            hMax.append( hDisp[peaks_h[i]])
            TimeMaxH.append( timeList[ peaks_h[i] ])
        elif hDisp[peaks_h[i]] < 0:
            hMax = [] #re-initialize list if a negative value occurs
            TimeMaxH = []

    
    for i in range(1, len(peaks_teta)): 
        if len(tetaMax) >= 5: 
            break #first if-statement to decide if there are enough data points
        elif tetaDisp[peaks_teta[i]] > 0:
            tetaMax.append( tetaDisp[peaks_teta[i]])
            TimeMaxTeta.append( timeList[ peaks_teta[i] ])
        elif tetaDisp[peaks_teta[i]] < 0:
            tetaMax = [] #re-initialize list if a negative value occurs
            TimeMaxTeta = []       


    zeta_h =    []          #liste d'amortissement en h
    zeta_teta = []          #liste d'amortissement en h
    delta_h = .0            #initialisation
    current_zeta_h = .0     #initialisation
    delta_teta = .0         #initialisation
    current_zeta_teta = .0  #initialisation
    
    """Again creating two separate loops"""
    for i in range(len(hMax)-1):
        delta_h = (np.log(hMax[i+1])-np.log(hMax[i]))
        current_zeta_h = delta_h / np.sqrt((2*PI)**2 + delta_h**2)
        zeta_h.append(current_zeta_h)
        
    for i in range(len(tetaMax)-1):
        delta_teta = (np.log(tetaMax[i+1])-np.log(tetaMax[i]))
        current_zeta_teta = delta_teta / np.sqrt((2*PI)**2 + delta_teta**2)
        zeta_teta.append(current_zeta_teta)
     
    avg_damping_h = .0
    if len(zeta_h) != 0:
        avg_damping_h = round(sum(zeta_h) / len(zeta_h), 3)
    else:
        print("damping_h list has no values")
    var_damping_h = float("{:.3e}".format(np.var(zeta_h)))
    avg_damping_teta = .0
    if len(zeta_teta) != 0:
        avg_damping_teta = round(sum(zeta_teta) / len(zeta_teta), 3)
    else:
        print("damping_teta list has no values")
    var_damping_teta = float("{:.3e}".format(np.var(zeta_teta)))
    
    data = [avg_damping_h, var_damping_h, avg_damping_teta, var_damping_teta]
    
    return data

def computeFrequencies (Nb, posH, posTeta, MaxTime):
    DispH2 = []
    DispTeta2 = []
    Ind1 = int(Nb / 2) - 1
    DispH2 = posH[Ind1:-1]
    DispTeta2 = posTeta[Ind1:-1]
    #
    FFT1 = abs(fftpack.fft(DispH2))
    FFTH = FFT1.tolist()
    FFT2 = abs(fftpack.fft(DispTeta2))
    FFTTeta = FFT2.tolist()
    #
    FFTH[0] = 0.
    FFTTeta[0] = 0.
    fH = 0.0
    fTeta = 0.0
    # max(FFTH[0:(np.size(FFTH) - 1) // 2 # Divide by 2 since FFT is symmetric 
    # So take also MaxT/2
    fH = (FFTH.index(max(FFTH[0:(np.size(FFTH) - 1) // 2]))) / (MaxTime / 2.)
    fTeta = (FFTTeta.index(max(FFTTeta[0:(np.size(FFTTeta) - 1) // 2]))) / (MaxTime / 2.)
    
    return fH, fTeta;