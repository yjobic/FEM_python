# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

from scipy.sparse.linalg import spsolve
import matplotlib.pyplot  as plt
import numpy as np
import time


#local import
from FEMlib.mesh import *
from FEMlib.assembleProblem import *
from FEMlib.quadrature import *
from FEMlib.basis import *
from FEMlib.plotSol import *
from FEMlib.errors import *

def g(x,y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)
def f(x,y):
    return g(x,y)*(2.*np.pi*np.pi)

def printTemps(temps):
    if temps>3600: 
        return '{:2.3f}'.format(temps/3600)+' h'
    elif temps>60: 
        return '{:2.3f}'.format(temps/60)+' min'
    else: 
        return '{:2.3f}'.format(temps)+' s'

outPlot=0
if  outPlot :
    get_ipython().run_line_magic('matplotlib', 'qt')
else :
    get_ipython().run_line_magic('matplotlib', 'inline')

start = time.time()

#load the mesh
mesh = Mesh()
startloc = time.time()
MeshFileName="Meshes/square_T2_1.msh"
mesh.GmshToMesh(MeshFileName)
endloc = time.time(); elapsed = endloc - startloc
print('Temps GmshToMesh : '+ printTemps(elapsed))

#On positionne les espaces d'approximations internes
mesh.updateBaseQuad(3)

startloc = time.time()
t=assembleMatrix(mesh, 1, 0)
endloc = time.time(); elapsed = endloc - startloc
print('Temps assembleMatrix : ' + printTemps(elapsed))

startloc = time.time()
B=assembleVecForceElem(mesh, f)
endloc = time.time(); elapsed = endloc - startloc
print('Temps assembleVecForceElem : ' + printTemps(elapsed))

Dirichlet(mesh,t,g,B)

A=coo_matrix(t.data).tocsr()
X = spsolve(A, B) # solve avec scipy
print(f'min(X) : {min(X)} - max(X) : {max(X)}')

end = time.time()
elapsed = end - start
print('Temps d\'exécution : ' + printTemps(elapsed))

SolEx=np.zeros(mesh.Npts)
for i in range(mesh.Npts):
    SolEx[i]=g(mesh.points[i].coord[0],mesh.points[i].coord[1])


# Compute errors
rmsE = rmsError(X, SolEx)
print(f'RMS error : {rmsE}')
errL2 = L2error(mesh,X,g)
print(f'L2 error : {errL2}')
print(f'Nombre total d elements : {mesh.getTotEle()}')
print(f'Nombre total de dof : {mesh.Npts}')

drawMesh=True
if (mesh.getTotEle()<500): 
    plotMesh(mesh)
else:
    drawMesh=False
plotSol(mesh,X,drawMesh)

