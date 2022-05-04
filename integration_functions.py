# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 18:20:11 2022

@author: Perseus - Raphaël AUBRY, Eliot MARTIN, Adrien SIMONNET
"""


def UpdateVP(MassTerm, CouplingTerm, DissTerm, StiffTerm, AeroTerm, XI):
    A = [0.0, 0.0]
    A[0] = XI[1]
    A[1] = 1/MassTerm * (AeroTerm - CouplingTerm - StiffTerm - DissTerm)
    # MassTerm = M*b (if ddh)
    # MassTerm = 1/r_alpha**2 (if ddteta)
    return A

def IntegrateRK4(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, Vec1, DT):
    X0 = [Vec1[0], Vec1[1]]
    X1 = [0.0, 0.0]
    X2 = [0.0, 0.0]
    X3 = [0.0, 0.0]
    A0 = [0.0, 0.0]
    A1 = [0.0, 0.0]
    A2 = [0.0, 0.0]
    A3 = [0.0, 0.0]

    A0 = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X0)

    X1[0] = X0[0] + A0[0] * DT / 2.0
    X1[1] = X0[1] + A0[1] * DT / 2.0

    A1 = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X1)
    X2[0] = X0[0] + A1[0] * DT / 2.0
    X2[1] = X0[1] + A2[1] * DT / 2.0

    A2 = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X2)
    X3[0] = X0[0] + A2[0] * DT / 1.0
    X3[1] = X0[1] + A2[1] * DT / 1.0

    A3 = UpdateVP(Mass, CouplingTerm, DissTerm, StiffTerm, AeroTerm, X3)

    Vec1[0] = Vec1[0] + (DT / 6.0) * (A0[0] + 2.0 * A1[0] + 2.0 * A2[0] + A3[0])
    Vec1[1] = Vec1[1] + (DT / 6.0) * (A0[1] + 2.0 * A1[1] + 2.0 * A2[1] + A3[1])
    Vec1[2] = (1. / 6.) * (A0[1] + 2.0 * A1[1] + 2.0 * A2[1] + A3[1])

    return Vec1
