# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

from FEMlib.Triplets import *
import numpy as np
from FEMlib.basis import *

def assembleMatrix(mesh, a, b):
    #tester sur idGroupEle est bien dans mesh
    
    #Matrice brut, non assemble
    t = Triplets()
    B = np.zeros(mesh.Npts)
    
    for idGroupEle in range(mesh.NbGroupEle):
    
        #On bosse sur chaque ele
        elelist=mesh.listElesType[idGroupEle]
        basis=mesh.base[idGroupEle]
        quad=mesh.quad[idGroupEle]
        MassElem= np.zeros((basis.Nphi,basis.Nphi))
        RigidElem= np.zeros((basis.Nphi,basis.Nphi))
        
        for ele in elelist:
            #print(ele)
            MassElem.fill(0)
            RigidElem.fill(0)
            #On construit les matrices elementaires dans l'espace de reference
            if b!= 0: 
                constructMassElem(basis, MassElem, ele, quad)
            if a!= 0: 
                constructRigidElem(basis, RigidElem, ele, quad)
            
            npEle=len(ele.p)
            for i in range(npEle):
                for j in range(npEle):
                    #coeff=a*RigidElem[i][j]*ele.vol+b*MassElem[i][j]
                    coeff=a*RigidElem[i][j]+b*MassElem[i][j]
                    t.append(ele.p[i].id, ele.p[j].id, coeff)
    
    return t

def assembleVecForceElem(mesh, fonc):
    #tester sur idGroupEle est bien dans mesh
    
    #Matrice brut, non assemble
    B = np.zeros(mesh.Npts)
    
    for idGroupEle in range(mesh.NbGroupEle):
    
        #On bosse sur chaque ele
        elelist=mesh.listElesType[idGroupEle]
        basis=mesh.base[idGroupEle]
        quad=mesh.quad[idGroupEle]
        
        for ele in elelist:
            #On construit les matrices elementaires dans l'espace de reference
            VecSolElem=constructVectSolElem(mesh, ele, basis, quad, fonc)
            
            npEle=len(ele.p)
            for i in range(npEle):
                B[ele.p[i].id] += VecSolElem[i]
    
    return B



def Dirichlet(mesh, t, fonc, B):
    #On se ballade sur tous les groupes de CL
    for gCL in range(len(mesh.CL)):
        #Dans un groupe de CL, on a les nums de dof (i.e. les num des pts)
        for pts in mesh.CL[gCL]:
            #On cherche dans la matrice les val a annuler
            #t.data[1][0] : indice i
            for dof in np.where(t.data[1][0]==pts):
                #On modifie la valeur de la contribution
                t.data[0][dof]=0
            t.append(int(pts),int(pts),1.)
            B[int(pts)]=fonc(mesh.points[pts].coord[0],mesh.points[pts].coord[1])











            