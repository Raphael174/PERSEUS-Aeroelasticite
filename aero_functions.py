# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 18:38:17 2022

@author: Perseus - Raphaël AUBRY, Eliot MARTIN, Adrien SIMONNET
"""

PI = 3.1415926535898

def computeCk(U_inf, b, omega_teta):
    """Changement de la fonction, usage de b au lieu de c"""

    k = omega_teta * b / U_inf
    Z0 = complex(1, (0.0455 / k))
    Z1 = complex(1, (0.3 / k))
    Z3 = 1 - (0.165 / Z0) - (0.335 / Z1)
    Ck = (Z3.real ** 2 + Z3.imag ** 2) ** 0.5
    return Ck

def computeQ(rho,U_inf):
    q_inf = 0.5*rho*U_inf**2
    return q_inf

def computeLift(omega_teta, rho, dCl, a, b, c, U_inf, Surf, VecH, VecTeta):
    
    #Implmentation of Claudia's work
    #q_inf = (rho*b*U_inf**2*Surf)/2
    #t_0 = b/U_inf
    
    
    #FA = 2*PI*q_inf*b * (VecTeta[1]*t_0 + VecH[2]*t_0**2/b -(c - a)*VecTeta[2]*t_0**2 )
    #FA += 2*dCl*computeCk(U_inf, b, omega_teta) * ( VecH[1]*t_0/b + VecTeta[0] + (0.5 - a)*t_0*VecTeta[1])
    # ------------------
    
    #Correction BE by Gourdain:

    CZ = -PI*( b/U_inf*VecTeta[1] + b/U_inf**2*VecH[2] - (a*b**2)/U_inf**2*VecTeta[2] )
    CZ += -dCl*computeCk(U_inf,c, omega_teta)*(VecTeta[0] + VecH[1]/U_inf + ( (0.5-a)*(b*VecTeta[1])/U_inf) )
    #------------------

    FA = computeQ(rho,U_inf)*Surf*(CZ) #compute force
    
    return (FA)

def computeMoment(omega_teta, rho, dCm, a, b, c, U_inf, Surf, VecH, VecTeta):
    #Implmentation of Claudia's work
    # q_inf = (rho*b*U_inf**2*Surf)/2
    # t_0 = b/U_inf


    # MA = 2*PI*q_inf*b**2 * (-(c-a)*t_0**2/b * VecH[2] + (0.5 - (c-a))*t_0 * VecTeta[1] + (0.165 + (c - a)**2)*t_0**2 * VecTeta[2])
    # MA += 2*dCm*b**2*q_inf*(0.5 + (c - a))*computeCk(U_inf, b, omega_teta) * (t_0/b + VecTeta[0] + (0.5 - (c - a))*t_0 * VecTeta[1])
    # ------------------

    #Correction BE by Gourdain:

    CM = PI/2*((a*b)/U_inf**2*VecH[2] - b**2/U_inf**2*(1./8.+a**2)*VecTeta[2] + (a-0.5)*b*VecTeta[1]/U_inf)
    CM += dCm*computeCk(U_inf,c, omega_teta)*(VecTeta[0] + VecH[1]/U_inf + ((0.5-a)*(b*VecTeta[1])/U_inf))

    MA = computeQ(rho,U_inf)*Surf*c*(CM); #compute moment

    return (MA)
