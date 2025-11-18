# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 16:27:37 2025

@author: yann
"""
import sys
#To install: conda install -c conda-forge gmsh python-gmsh
import gmsh
#Doc API gmsh : https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_12_1/api/gmsh.py

import matplotlib.pyplot  as plt
from matplotlib.tri import Triangulation

import os, sys
import numpy as np

sys.path.insert(0, os.path.realpath('../'))

from FEMlib.mesh import *
from FEMlib.plotSol import *


np.random.seed(10)

OpenFile = "square_T1_2.msh"
SavedFile = "square_T1_2_displaced.msh"
Amplitude = 0.07

gmsh.initialize()
gmsh.open(OpenFile)

# Tous les noeuds
nodeTags, nodeCoords, paramCoords = gmsh.model.mesh.getNodes()

# Obtention des noeuds frontières (par exemple pour 2D, dim=1)
edge_nodes = set()
for dim, tag in gmsh.model.getPhysicalGroups():
    if dim == 1:
        tags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
        edge_nodes.update(tags)

# Déplacement random des noeuds internes
for i, tag in enumerate(nodeTags):
    if tag not in edge_nodes:
        x = nodeCoords[3*i]
        y = nodeCoords[3*i + 1]
        z = nodeCoords[3*i + 2]
        dx, dy = np.random.uniform(-Amplitude, Amplitude, 2)
        newPos = [x + dx, y + dy, z]
        pCoord = paramCoords[2*i:2*i+2] if len(paramCoords) > 0 else []
        gmsh.model.mesh.setNode(tag, newPos, pCoord)

gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.write(SavedFile)
gmsh.finalize()

mesh = Mesh()
mesh.GmshToMesh(SavedFile)
plotMesh(mesh)
