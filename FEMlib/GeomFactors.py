# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

import numpy as np


def computeJac(ele):
    Jac=np.array([[np.double(0),np.double(0)],[np.double(0),np.double(0)]])
    for i in range(2):
        for j in range(2):
            #ici ele triangulaire, forme analytique
            Jac[i][j]=ele.p[j+1].coord[i]-ele.p[0].coord[i]
    return Jac

def computeJacGeneral(base, quad, ele, Jac, xi, eta):
    """
    On calcule la jacobienne au point de quadrature k
        -------
    Jac : Array
    """
    #i -> (dxi,deta)
    for i in range(2): #2 : c'est la dim physique 2D ou 3D
        #j -> (x,y)
        for j in range(2): #2 : c'est la dim physique 2D ou 3D
            for N in range(base.Nphi):
                Jac[i][j]+=base.DeriveesFctForm[N][i](xi,eta)*ele.p[N].coord[j]


def computeCoordPhyFromRef(ele, base, quad):
    x=np.zeros(quad.Npts)
    y=np.zeros(quad.Npts)
    #On cherche la position des points de gauss dans l'espace physique
    for k in range(quad.Npts):
        xi= quad.coords[k][0]
        eta=quad.coords[k][1]
        for i in range(base.Nphi):
            x[k]+=base.FctForm[i](xi,eta)*ele.p[i].coord[0]
            y[k]+=base.FctForm[i](xi,eta)*ele.p[i].coord[1]
    return x,y

            
def computeDetJ(Jac):
    return np.linalg.det(Jac)

def computeInvJ(Jac):
    return np.linalg.inv(Jac)
    
  
