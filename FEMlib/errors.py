# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

from FEMlib.mesh import *


def rmsError(Sol, SolEx):
    relative_error=np.sum(np.power(Sol - SolEx,2))
    return np.sqrt(relative_error/len(Sol))

def interpoleSolAtGaussPoints(ele,quad,base,Sol):
    SolInterpole=np.zeros(quad.Npts)
    for k in range(quad.Npts):
        xi= quad.coords[k][0]
        eta=quad.coords[k][1]
        for i in range(base.Nphi):
            SolInterpole[k]+=base.FctForm[i](xi,eta)*Sol[ele.p[i].id]
    return SolInterpole

def errorL2Elem(ele,quad,base,Sol,FuncSolEx):
    errorElem=0
    Jac=np.array([[np.double(0),np.double(0)],[np.double(0),np.double(0)]])
    SolInterpole=interpoleSolAtGaussPoints(ele,quad,base,Sol)
    x,y=computeCoordPhyFromRef(ele, base, quad)      
    for k in range(quad.Npts):
        xi =quad.coords[k][0]
        eta=quad.coords[k][1]
        #Jacobian is not constant, we have to evalute it at gauss point
        computeJacGeneral(base, quad, ele, Jac, xi, eta)
        coeff = computeDetJ(Jac)*quad.wp[k]
        errorElem+=coeff * pow(SolInterpole[k] - FuncSolEx(x[k],y[k]),2)
    return errorElem

def L2error(mesh, Sol, FuncSolEx):
    errL2=0
    for idGroupEle in range(mesh.NbGroupEle):
        elelist=mesh.listElesType[idGroupEle]
        basis=mesh.base[idGroupEle]
        quad=mesh.quad[idGroupEle]
        for ele in elelist:
            errL2+=errorL2Elem(ele,quad,basis,Sol,FuncSolEx)
    return np.sqrt(errL2)