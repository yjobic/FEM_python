# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

from scipy.sparse.linalg import spsolve
import matplotlib.pyplot  as plt
import numpy as np
import time
import sys


#local import
from FEMlib.mesh import *
from FEMlib.assembleProblem import *
from FEMlib.quadrature import *
from FEMlib.basis import *
from FEMlib.plotSol import *
from FEMlib.errors import *

def g(x,y):
    #
    #[A REMPLIR] : solution exacte, utilisée pour les conditions de Dirichlet
    #
    return None  # à remplacer

def f(x,y):
    #
    #[A REMPLIR] : terme source (second membre)
    #
    return None  # à remplacer

# Vérification que g et f ont bien été complétées avant de lancer le calcul
fonctions_a_remplir = [nom for nom, fct in (('g', g), ('f', f))
                       if fct(0.5, 0.5) is None]
if fonctions_a_remplir:
    print('\n' + '*'*70)
    print(' ERREUR : fonction(s) à compléter : ' + ', '.join(fonctions_a_remplir))
    print(' Remplissez les parties [A REMPLIR] au début de solve_edp.py :')
    print('   - g(x,y) : la solution exacte (conditions de Dirichlet)')
    print('   - f(x,y) : le terme source (second membre)')
    print(' puis relancez le programme.')
    print('*'*70 + '\n')
    sys.exit(1)

def printTemps(temps):
    if temps>3600: 
        return '{:2.3f}'.format(temps/3600)+' h'
    elif temps>60: 
        return '{:2.3f}'.format(temps/60)+' min'
    else: 
        return '{:2.3f}'.format(temps)+' s'

outPlot=0
try: # Vérification sécurisée de l'environnement IPython
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        if outPlot:
            ipython.run_line_magic('matplotlib', 'qt')
        else:
            ipython.run_line_magic('matplotlib', 'inline')
except ImportError:
    pass  # IPython n'est pas installé ou est indisponible
    
start = time.time()

#load the mesh
mesh = Mesh()
startloc = time.time()
MeshFileName="Meshes/square_Q_1_o2.msh"
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

