# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

import sys

#local import
from FEMlib.Triplets import *
from FEMlib.GeomFactors import *


def constructMassElem(base, matrix, ele, quad):
    Jac=np.array([[np.double(0),np.double(0)],[np.double(0),np.double(0)]])
    for k in range(quad.Npts):
        Jac.fill(0)
        xi=quad.coords[k][0]
        eta=quad.coords[k][1]
        computeJacGeneral(base, quad, ele, Jac, xi, eta)
        coeff = computeDetJ(Jac)*quad.wp[k]
        for i in range(base.Nphi):
            for j in range(base.Nphi):
                matrix[i,j] += coeff \
                              *base.FctForm[i](xi,eta) \
                              *base.FctForm[j](xi,eta)

def constructRigidElem(base, matrix, ele, quad):
    Jac=np.array([[np.double(0),np.double(0)],[np.double(0),np.double(0)]])
    for k in range(quad.Npts):
        Jac.fill(0)
        xi=quad.coords[k][0]
        eta=quad.coords[k][1]
        computeJacGeneral(base, quad, ele, Jac, xi, eta)
        Bp = computeInvJ(Jac)
        coeff = computeDetJ(Jac)*quad.wp[k]
        #coeff = quad.wp[k]
        for i in range(base.Nphi):
            phi1=np.array([base.DeriveesFctForm[i][0](xi,eta),base.DeriveesFctForm[i][1](xi,eta)])
            gphii=np.matmul(Bp,phi1)
            for j in range(base.Nphi):
                phi2=np.array([base.DeriveesFctForm[j][0](xi,eta),base.DeriveesFctForm[j][1](xi,eta)])
                gphij=np.matmul(Bp,phi2)
                val=coeff*np.matmul(gphii,gphij)
                matrix[i,j]+=val

def constructVectSolElem(mesh, ele, base, quad, fonc):
    Jac=np.array([[np.double(0),np.double(0)],[np.double(0),np.double(0)]])
    vect = np.zeros(len(ele.p))                
    x,y=computeCoordPhyFromRef(ele, base, quad)      
    for k in range(quad.Npts):
        Jac.fill(0)
        xi =quad.coords[k][0]
        eta=quad.coords[k][1]
        computeJacGeneral(base, quad, ele, Jac, xi, eta)
        coeff = computeDetJ(Jac)*quad.wp[k]
        #coeff = quad.wp[k]
        for i in range(base.Nphi):
            phi1=base.FctForm[i](xi,eta) 
            vect[i]+=coeff * phi1 * fonc(x[k],y[k])
    return vect


class Basis:
    def __init__():
        self.dim=0
        self.Nphi=0  #On a Nphi fonction de forme pour ndof
    
class LagrangeT1(Basis):
    name='Base T1 lagrange \n'
    #Datt thouzot p 108 (2.3.2)

    def __init__(self,dim,ndof):
        self.dim=dim
        self.Nphi=ndof #On a Nphi fonction de forme pour ndof
        assert dim==2, "Attention on ne gere que dim=2"
        
    def Phi1(xi,eta):
        return 1-eta-xi
    def Phi2(xi,eta):
        return xi
    def Phi3(xi,eta):
        return eta
    FctForm=[Phi1,Phi2,Phi3]
    
    #On passe aux derivées des fonctions de formes dans l'espace de reference
    #On passe aux derivées des fonctions de formes
    def dPhi1dxi(xi,eta):
        return -1
    def dPhi1deta(xi,eta):
        return -1
    def dPhi2dxi(xi,eta):
        return 1
    def dPhi2deta(xi,eta):
        return 0
    def dPhi3dxi(xi,eta):
        return 0
    def dPhi3deta(xi,eta):
        return 1
    DeriveesFctForm=[[dPhi1dxi,dPhi1deta],[dPhi2dxi,dPhi2deta],[dPhi3dxi,dPhi3deta]]                     


class LagrangeT2(Basis):
    name='Base T2 lagrange\n'
    #Datt thouzot p 110 (2.3.3.1)
    
    def __init__(self,dim,ndof):
        self.dim=dim
        self.Nphi=ndof #On a Nphi fonction de forme pour ndof
        assert dim==2, "Attention on ne gere que dim=2"
        
    def Phi1(xi,eta):
        buf=1-xi-eta
        return -buf*(1-2*buf)
    def Phi2(xi,eta):
        buf=1-xi-eta
        return 4*xi*buf
    def Phi3(xi,eta):
        return -xi*(1-2*xi)
    def Phi4(xi,eta):
        return 4*xi*eta
    def Phi5(xi,eta):
        return -eta*(1-2*eta)
    def Phi6(xi,eta):
        buf=1-xi-eta
        return 4*eta*buf
    FctForm=[Phi1,Phi2,Phi3,Phi4,Phi5,Phi6]
        
    #On passe aux derivées des fonctions de formes
    def dPhi1dxi(xi,eta):
        buf=1-xi-eta
        return 1-4*buf
    def dPhi1deta(xi,eta):
        buf=1-xi-eta
        return 1-4*buf
    def dPhi2dxi(xi,eta):
        buf=1-xi-eta
        return 4*(buf-xi)
    def dPhi2deta(xi,eta):
        return -4*xi
    def dPhi3dxi(xi,eta):
        return -1+4*xi
    def dPhi3deta(xi,eta):
        return 0
    def dPhi4dxi(xi,eta):
        return 4*eta
    def dPhi4deta(xi,eta):
        return 4*xi
    def dPhi5dxi(xi,eta):
        return 0
    def dPhi5deta(xi,eta):
        return -1+4*eta
    def dPhi6dxi(xi,eta):
        return -4*eta
    def dPhi6deta(xi,eta):
        buf=1-xi-eta
        return 4*(buf-eta)
    DeriveesFctForm=[[dPhi1dxi,dPhi1deta],[dPhi2dxi,dPhi2deta],[dPhi3dxi,dPhi3deta], \
                     [dPhi4dxi,dPhi4deta],[dPhi5dxi,dPhi5deta],[dPhi6dxi,dPhi6deta]]



class LagrangeQ1(Basis):
    name='Base Q1 lagrange\n'
    #Datt thouzot p 121 (2.4.2)
    
    def __init__(self,dim,ndof):
        self.dim=dim
        self.Nphi=ndof #On a Nphi fonction de forme pour ndof
        assert dim==2, "Attention on ne gere que dim=2"
        
    def Phi1(xi,eta):
        return 0.25*(1-xi)*(1-eta)
    def Phi2(xi,eta):
        return 0.25*(1+xi)*(1-eta)
    def Phi3(xi,eta):
        return 0.25*(1+xi)*(1+eta)
    def Phi4(xi,eta):
        return 0.25*(1-xi)*(1+eta)
    FctForm=[Phi1,Phi2,Phi3,Phi4]
        
    #On passe aux derivées des fonctions de formes
    def dPhi1dxi(xi,eta):
        return 0.25*(eta-1)
    def dPhi1deta(xi,eta):
        return 0.25*(xi-1)
    def dPhi2dxi(xi,eta):
        return 0.25*(1-eta)
    def dPhi2deta(xi,eta):
        return 0.25*(-1-xi)
    def dPhi3dxi(xi,eta):
        return 0.25*(1+eta)
    def dPhi3deta(xi,eta):
        return 0.25*(1+xi)
    def dPhi4dxi(xi,eta):
        return 0.25*(-1-eta)
    def dPhi4deta(xi,eta):
        return 0.25*(1-xi)
    DeriveesFctForm=[[dPhi1dxi,dPhi1deta],[dPhi2dxi,dPhi2deta],[dPhi3dxi,dPhi3deta], \
                     [dPhi4dxi,dPhi4deta]]
        
class LagrangeQ2(Basis):
    name='Base Q2 lagrange\n'
    #Datt thouzot p 122 (2.4.3.1)
    
    def __init__(self,dim,ndof):
        self.dim=dim
        self.Nphi=ndof #On a Nphi fonction de forme pour ndof
        assert dim==2, "Attention on ne gere que dim=2"
        
    def Phi1(xi,eta):
        return 0.25*(1-xi)*(1-eta)*xi*eta
    def Phi2(xi,eta):
        return -0.5*(1-xi*xi)*(1-eta)*eta
    def Phi3(xi,eta):
        return -0.25*(1+xi)*(1-eta)*xi*eta
    def Phi4(xi,eta):
        return 0.5*(1+xi)*(1-eta*eta)*xi
    def Phi5(xi,eta):
        return 0.25*(1+xi)*(1+eta)*xi*eta
    def Phi6(xi,eta):
        return 0.5*(1-xi*xi)*(1+eta)*eta
    def Phi7(xi,eta):
        return -0.25*(1-xi)*(1+eta)*xi*eta
    def Phi8(xi,eta):
        return -0.5*(1-xi)*(1-eta*eta)*xi
    def Phi9(xi,eta):
        return (1-xi*xi)*(1-eta*eta)
    FctForm=[Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7,Phi8,Phi9]
        
    #On passe aux derivées des fonctions de formes
    def dPhi1dxi(xi,eta):
        return 0.25*(1-2*xi)*(1-eta)*eta
    def dPhi1deta(xi,eta):
        return 0.25*(1-xi)*(1-2*eta)*xi
    def dPhi2dxi(xi,eta):
        return (1-eta)*xi*eta
    def dPhi2deta(xi,eta):
        return -0.5*(1-xi*xi)*(1-2*eta)
    def dPhi3dxi(xi,eta):
        return -0.25*(1+2*xi)*(1-eta)*eta
    def dPhi3deta(xi,eta):
        return -0.25*(1+xi)*(1-2*eta)*xi
    def dPhi4dxi(xi,eta):
        return 0.5*(1+2*xi)*(1-eta*eta)
    def dPhi4deta(xi,eta):
        return -(1+xi)*xi*eta
    def dPhi5dxi(xi,eta):
        return 0.25*(1+2*xi)*(1+eta)*eta
    def dPhi5deta(xi,eta):
        return 0.25*(1+xi)*(1+2*eta)*xi
    def dPhi6dxi(xi,eta):
        return -(1+eta)*xi*eta
    def dPhi6deta(xi,eta):
        return 0.5*(1-xi*xi)*(1+2*eta)
    def dPhi7dxi(xi,eta):
        return -0.25*(1-2*xi)*(1+eta)*eta
    def dPhi7deta(xi,eta):
        return -0.25*(1-xi)*(1+2*eta)*xi
    def dPhi8dxi(xi,eta):
        return -0.5*(1-2*xi)*(1-eta*eta)
    def dPhi8deta(xi,eta):
        return (1-xi)*xi*eta
    def dPhi9dxi(xi,eta):
        return -2*(1-eta*eta)*xi
    def dPhi9deta(xi,eta):
        return -2*(1-xi*xi)*eta
    DeriveesFctForm=[[dPhi1dxi,dPhi1deta],[dPhi2dxi,dPhi2deta],[dPhi3dxi,dPhi3deta], \
                     [dPhi4dxi,dPhi4deta],[dPhi5dxi,dPhi5deta],[dPhi6dxi,dPhi6deta], \
                     [dPhi7dxi,dPhi7deta],[dPhi8dxi,dPhi8deta],[dPhi9dxi,dPhi9deta]]


















