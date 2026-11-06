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

ordre = 1
OpenFile = "square_T_1_o"+str(ordre)+".msh"
SavedFile = "square_T1_"+str(ordre)+"_displaced.msh"
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

# Coordonnées de chaque noeud, indexées par leur tag
coords = {tag: np.array(nodeCoords[3*i:3*i+3]) for i, tag in enumerate(nodeTags)}

# Éléments surfaciques (dim=2) : on distingue les sommets (noeuds "primaires")
# des noeuds d'ordre supérieur (milieux d'arêtes, centre des quadrangles)
elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2)
elements = []      # liste de (nombre de sommets, noeuds de l'élément)
vertex_nodes = set()
for etype, enodes in zip(elemTypes, elemNodeTags):
    props = gmsh.model.mesh.getElementProperties(etype)
    nNodes, nVertices = props[3], props[5]
    for e in np.reshape(enodes, (-1, nNodes)):
        elements.append((nVertices, e))
        vertex_nodes.update(e[:nVertices])

# 1) Déplacement random des sommets internes uniquement
for tag in sorted(vertex_nodes):
    if tag not in edge_nodes:
        dx, dy = np.random.uniform(-Amplitude, Amplitude, 2)
        coords[tag] = coords[tag] + np.array([dx, dy, 0.])

# 2) Repositionnement des noeuds d'ordre 2 à partir des sommets déplacés :
#    noeud milieu = centre du brin, noeud central (quad 9) = barycentre.
#    Numérotation gmsh : sommets, puis milieux des arêtes (s0s1, s1s2, ...),
#    puis éventuellement le noeud central.
for nV, e in elements:
    for k in range(nV, len(e)):
        if e[k] in edge_nodes:
            continue
        if k < 2*nV:  # milieu de l'arête (s_j, s_j+1)
            j = k - nV
            coords[e[k]] = 0.5 * (coords[e[j]] + coords[e[(j + 1) % nV]])
        else:         # noeud central
            coords[e[k]] = np.mean([coords[v] for v in e[:nV]], axis=0)

# Mise à jour des noeuds internes dans gmsh
for i, tag in enumerate(nodeTags):
    if tag not in edge_nodes:
        pCoord = paramCoords[2*i:2*i+2] if len(paramCoords) > 0 else []
        gmsh.model.mesh.setNode(tag, list(coords[tag]), pCoord)

gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.write(SavedFile)
gmsh.finalize()

# Comparaison maillage d'origine / maillage déformé
mesh_orig = Mesh()
mesh_orig.GmshToMesh(OpenFile)
mesh = Mesh()
mesh.GmshToMesh(SavedFile)

x0 = np.array([p.coord[0] for p in mesh_orig.points])
y0 = np.array([p.coord[1] for p in mesh_orig.points])
x1 = np.array([p.coord[0] for p in mesh.points])
y1 = np.array([p.coord[1] for p in mesh.points])

plt.figure(figsize=(8, 8))
# Maillage d'origine en gris pointillé, maillage déformé en noir
plotMesh(mesh_orig, show=False, edgecolor='grey', linestyle='--', linewidth=0.8)
plotMesh(mesh, show=False)
plt.plot(x0, y0, 'o', color='grey', ms=4)
plt.plot(x1, y1, 'ro', ms=4)
# Flèches : déplacement de chaque noeud
plt.quiver(x0, y0, x1 - x0, y1 - y0, angles='xy', scale_units='xy', scale=1,
           color='tab:blue', width=0.004)
# Légende (entrées "fantômes" pour les styles de traits)
plt.plot([], [], '--', color='grey', label="Maillage d'origine")
plt.plot([], [], '-', color='black', label="Maillage déformé")
plt.plot([], [], '-', color='tab:blue', label="Déplacement des noeuds")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)
plt.title(f"Noeuds internes déplacés (amplitude = {Amplitude})")
plt.tight_layout()
plt.show()
