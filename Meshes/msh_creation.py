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


def createMesh_triangle(filename, basesize, elementOrder):

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
    

def createMesh_square(filename, basesize, elementOrder):

    gmsh.initialize()
    gmsh.model.add("square_Q1")
    
    # Paramètres
    boxdim = 1.0
    gridsize = 1
    numTransfinitSize = int(1/basesize)+1
        
    # Géométrie "classique" (non-OCC) pour coller au .geo
    gmsh.model.geo.addPoint(0.0,      0.0,      0.0, gridsize, 1)
    gmsh.model.geo.addPoint(boxdim,   0.0,      0.0, gridsize, 2)
    gmsh.model.geo.addPoint(boxdim,   boxdim,   0.0, gridsize, 3)
    gmsh.model.geo.addPoint(0.0,      boxdim,   0.0, gridsize, 4)
    
    gmsh.model.geo.addLine(1, 2, 7)
    gmsh.model.geo.addLine(2, 3, 8)
    gmsh.model.geo.addLine(3, 4, 9)
    gmsh.model.geo.addLine(4, 1, 10)
    
    gmsh.model.geo.addCurveLoop([7, 8, 9, 10], 14)
    gmsh.model.geo.addPlaneSurface([14], 16)
    
    # Transfinite + recombine
    for l in [7, 8, 9, 10]:
        gmsh.model.geo.mesh.setTransfiniteCurve(l, numTransfinitSize)
    
    gmsh.model.geo.mesh.setTransfiniteSurface(16)
    gmsh.model.geo.mesh.setRecombine(2, 16)  # dim=2, tag=16
    
    # Synchronisation géométrie -> modèle
    gmsh.model.geo.synchronize()    

    # Groupes physiques
    # Limites 1D (tags 7–10)
    gmsh.model.addPhysicalGroup(1, [7, 8], 1)
    gmsh.model.setPhysicalName(1, 1, "Boundary 1")
    
    gmsh.model.addPhysicalGroup(1, [9, 10], 2)
    gmsh.model.setPhysicalName(1, 2, "Boundary 2")
    
    # Surface 2D
    gmsh.model.addPhysicalGroup(2, [16], 3)
    gmsh.model.setPhysicalName(2, 3, "Surface Rect")
    
    gmsh.option.setNumber("Mesh.Algorithm", 5)
        
    # Génération et sauvegarde
    gmsh.model.mesh.generate(2)
    
    gmsh.model.mesh.setOrder(elementOrder)

    gmsh.model.Format="msh41"
    gmsh.write(filename+".msh")
    gmsh.finalize()


order = 2
basesize = 0.5/2/2/2/2
MeshFileName="square_T_5_o"+str(order)
createMesh_triangle(MeshFileName,basesize,order)

mesh = Mesh()
mesh.GmshToMesh(MeshFileName+".msh")

plotMesh(mesh)


MeshFileName="square_Q_5_o"+str(order)
createMesh_square(MeshFileName,basesize,order)

mesh = Mesh()
mesh.GmshToMesh(MeshFileName+".msh")

plotMesh(mesh)
