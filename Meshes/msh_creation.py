# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

import sys
#To install: conda install -c conda-forge gmsh python-gmsh
import gmsh
#Doc API gmsh : https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_12_1/api/gmsh.py

import matplotlib.pyplot  as plt
from matplotlib.tri import Triangulation

import os, sys

sys.path.insert(0, os.path.realpath('../'))

from FEMlib.mesh import *
from FEMlib.plotSol import *


def createMesh(filename, basesize, elementOrder):

    gmsh.initialize(sys.argv)
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", basesize);
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", basesize);
    # Model
    model = gmsh.model
    model.add("Square")
    # Rectangle of (elementary) tag 1
    factory = model.occ
    factory.addRectangle(0,0,0, 1, 1, 1)
    # Sync
    factory.synchronize()
    # Physical groups
    gmsh.model.addPhysicalGroup(1, [1], 1)
    gmsh.model.addPhysicalGroup(1, [2,3,4], 2)
    gmsh.model.addPhysicalGroup(2, [1], 10)
    # Mesh (2D)
    model.mesh.generate(2)
    #order of the mesh, after generate
    gmsh.model.mesh.setOrder(elementOrder)
    #gmsh.model.mesh.HighOrderOptimize(2)

    #Save mesh
    model.Format="msh41"
    gmsh.write(filename+".msh")    
    gmsh.finalize()
    

MeshFileName="square_tri_fin_fin_fin_o1"
createMesh(MeshFileName,0.5/2/2/2,1)

mesh = Mesh()
mesh.GmshToMesh(MeshFileName+".msh")

plotMesh(mesh)
