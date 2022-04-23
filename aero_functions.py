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


def computeLift(omega_teta, rho, a, b, c, U_inf, Surf, half_c, VecH, VecTeta):
    
    dCl = 2*PI
    q_inf = (rho*half_c*U_inf**2*Surf)/2
    t_0 = b/U_inf
    
    
    FA = 2*PI*q_inf*b * (VecTeta[1]*t_0 + VecH[2]*t_0**2/b -(c - a)*VecTeta[2]*t_0**2 )
    FA += 2*dCl*computeCk(U_inf, b, omega_teta) * ( VecH[1]*t_0/b + VecTeta[0] + (0.5 - a)*t_0*VecTeta[1])
    # ------------------

    return FA

def computeMoment(omega_teta, rho, a, b, c, U_inf, Surf, half_c, VecH, VecTeta):

    dCm = 2*PI
    q_inf = (rho*half_c*U_inf**2*Surf)/2
    t_0 = b/U_inf


    MA = 2*PI*q_inf*b**2 * (-(c-a)*t_0**2/b * VecH[2] + (0.5 - (c-a))*t_0 * VecTeta[1] + (0.165 + (c - a)**2)*t_0**2 * VecTeta[2])
    MA += 2*dCm*b**2*q_inf*(0.5 + (c - a))*computeCk(U_inf, b, omega_teta) * (t_0/b + VecTeta[0] + (0.5 - (c - a))*t_0 * VecTeta[1])
    # ------------------

    return MA