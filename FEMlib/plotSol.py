# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

import matplotlib.tri as tri
import matplotlib.pyplot  as plt
import numpy as np
import sys

#local import
from FEMlib.mesh import *

#inspired by https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib

#***************************************************************
#  Partie concernant l'affichage de la solution
#***************************************************************

def interpolQuadCentre(base, quad, sol):
    #Centre du rectangle est en (xi,eta)=(0,0)
    #car rectangle va de [-1,1]^2
    xi=0; eta=0;
    val=0
    for i in range(len(quad)):
        val+=base.FctForm[i](xi,eta)*sol[quad[i]]
    return val

def ajoutBarycentre(x,y,quad):
    xb=0; yb=0;
    pos=len(x)
    for i in range(4):
        xb+=x[quad[i]]
        yb+=y[quad[i]]
    x.append(xb/4)
    y.append(yb/4)
    return pos


# converts quad elements into tri elements
def quads_to_tris(nodes_x, nodes_y, base, sol, quads):
    if base is not None:
        id(sol)
        newsol = np.zeros(len(sol)+len(quads))
        newsol[0:len(sol)]=sol
    tris = [[None for j in range(3)] for i in range(4*len(quads))]
    for i in range(len(quads)):
        j = 4*i
        n0 = quads[i][0]
        n1 = quads[i][1]
        n2 = quads[i][2]
        n3 = quads[i][3]
        n4 = ajoutBarycentre(nodes_x, nodes_y, quads[i])
        if base is not None:
            val=interpolQuadCentre(base, quads[i], sol)
            newsol[len(sol)+i] = val
        
        tris[j][0] = n0
        tris[j][1] = n4
        tris[j][2] = n3
        
        tris[j + 1][0] = n0
        tris[j + 1][1] = n1
        tris[j + 1][2] = n4
        
        tris[j + 2][0] = n1
        tris[j + 2][1] = n2
        tris[j + 2][2] = n4
        
        tris[j + 3][0] = n4
        tris[j + 3][1] = n2
        tris[j + 3][2] = n3
    if base is not None:
        return tris, newsol
    else:
        return tris

def triOrdre2_to_TrisOrdre1(tris_o2):
    tris_o1 = [[None for j in range(3)] for i in range(4*len(tris_o2))]
    big_tris_o1 = [[None for j in range(3)] for i in range(len(tris_o2))]
    for i in range(len(tris_o2)):
        j=4*i
        n0 = tris_o2[i][0]
        n1 = tris_o2[i][1]
        n2 = tris_o2[i][2]
        n3 = tris_o2[i][3]
        n4 = tris_o2[i][4]
        n5 = tris_o2[i][5]
        
        tris_o1[j][0] = n0
        tris_o1[j][1] = n1
        tris_o1[j][2] = n5

        tris_o1[j+1][0] = n1
        tris_o1[j+1][1] = n2
        tris_o1[j+1][2] = n3

        tris_o1[j+2][0] = n5
        tris_o1[j+2][1] = n1
        tris_o1[j+2][2] = n3

        tris_o1[j+3][0] = n5
        tris_o1[j+3][1] = n3
        tris_o1[j+3][2] = n4
        
        big_tris_o1[i][0] = n0
        big_tris_o1[i][1] = n2
        big_tris_o1[i][2] = n4
    return tris_o1,big_tris_o1
        
def quadsOrdre2_to_QuadsOrdre1(quads_o2):
    quads_o1 = [[None for j in range(4)] for i in range(4*len(quads_o2))]
    big_quads_o1 = [[None for j in range(4)] for i in range(len(quads_o2))]
    for i in range(len(quads_o2)):
        j=4*i
        n0 = quads_o2[i][0]
        n1 = quads_o2[i][1]
        n2 = quads_o2[i][2]
        n3 = quads_o2[i][3]
        n4 = quads_o2[i][4]
        n5 = quads_o2[i][5]
        n6 = quads_o2[i][6]
        n7 = quads_o2[i][7]
        n8 = quads_o2[i][8]
       
        quads_o1[j][0] = n0
        quads_o1[j][1] = n1
        quads_o1[j][2] = n8
        quads_o1[j][3] = n7

        quads_o1[j+1][0] = n1
        quads_o1[j+1][1] = n2
        quads_o1[j+1][2] = n3
        quads_o1[j+1][3] = n8

        quads_o1[j+2][0] = n8
        quads_o1[j+2][1] = n3
        quads_o1[j+2][2] = n4
        quads_o1[j+2][3] = n5

        quads_o1[j+3][0] = n7
        quads_o1[j+3][1] = n8
        quads_o1[j+3][2] = n5
        quads_o1[j+3][3] = n6
        
        big_quads_o1[i][0] = n0
        big_quads_o1[i][1] = n2
        big_quads_o1[i][2] = n4
        big_quads_o1[i][3] = n6
    return quads_o1,big_quads_o1


        
# plots a finite element mesh
def plot_fem_mesh(nodes_x, nodes_y, elements):
    for element in elements:
        x = [nodes_x[element[i]] for i in range(len(element))]
        y = [nodes_y[element[i]] for i in range(len(element))]
        plt.fill(x, y, edgecolor='black', fill=False)


def plotMesh(mesh):

    #Coordonnees des points
    x= [pt.coord[0] for pt in mesh.points]
    y= [pt.coord[1] for pt in mesh.points]

    #On construit la liste des noeuds composant les elements tri et quad
    triangle=[]
    quadrangle=[]
    triangle_o1=[]
    quadrangle_o1=[]
    big_triangle_o1=[]
    big_quadrangle_o1=[]
    for GrpEle in range(len(mesh.listElesType)):
        elelist=mesh.listElesType[GrpEle]
        if mesh.getNumNodePerEle(GrpEle) == 3: #C'est des triangles d'ordre 1
            for ltri in elelist:
                triangle.append([ p.id for p in ltri.p])
        elif mesh.getNumNodePerEle(GrpEle) == 4:
            for quad in elelist:
                quadrangle.append([ p.id for p in quad.p])
        elif mesh.getNumNodePerEle(GrpEle) == 6:
            triangle_o2=[]
            for ltri in elelist:
                triangle_o2.append([ p.id for p in ltri.p])
            triangle_o1,big_triangle_o1=triOrdre2_to_TrisOrdre1(triangle_o2)
        elif mesh.getNumNodePerEle(GrpEle) == 9:
            quads_o2=[]
            for quad in elelist:
                quads_o2.append([ p.id for p in quad.p])
            quadrangle_o1,big_quadrangle_o1=quadsOrdre2_to_QuadsOrdre1(quads_o2)
        else:
            print('\nNombre de noeud par element non gere : ',mesh.getNumNodePerEle(GrpEle))
            sys.exit(1)
            
    elements = triangle + quadrangle + big_triangle_o1 + big_quadrangle_o1

    # convert all elements into triangles
    elements_all_tris = triangle + quads_to_tris(x, y, None, None, quadrangle) \
                      + triangle_o1 + quads_to_tris(x,y, None, None, quadrangle_o1)

    # create an unstructured triangular grid instance
    triangulation = tri.Triangulation(x, y, elements_all_tris)

    # plot the finite element mesh
    plot_fem_mesh(x, y, elements)

    # show
    plt.axis('equal')
    plt.show()

def plotSol(mesh, X, printMesh):

    #Coordonnees des points
    x= [pt.coord[0] for pt in mesh.points]
    y= [pt.coord[1] for pt in mesh.points]

    #On construit la liste des noeuds composant les elements tri et quad
    triangle=[]
    quadrangle=[]
    triangle_o1=[]
    quadrangle_o1=[]
    big_triangle_o1=[]
    big_quadrangle_o1=[]
    tri1=[]
    for GrpEle in range(len(mesh.listElesType)):
        NumNodePerEle=mesh.getNumNodePerEle(GrpEle)
        if NumNodePerEle == 3 or NumNodePerEle == 6 :
            base=LagrangeT1(2,3)
        elif NumNodePerEle == 4 or NumNodePerEle == 9:
            base=LagrangeQ1(2,4)
        elelist=mesh.listElesType[GrpEle]
        if mesh.getNumNodePerEle(GrpEle) == 3: #C'est des triangles d'ordre 1
            for ltri in elelist:
                triangle.append([ p.id for p in ltri.p])
        elif mesh.getNumNodePerEle(GrpEle) == 4:
            for quad in elelist:
                quadrangle.append([ p.id for p in quad.p])
            tri1, X = quads_to_tris(x,y,base,X,quadrangle)
        elif mesh.getNumNodePerEle(GrpEle) == 6:
            triangle_o2=[]
            for ltri in elelist:
                triangle_o2.append([ p.id for p in ltri.p])
            triangle_o1,big_triangle_o1=triOrdre2_to_TrisOrdre1(triangle_o2)
        elif mesh.getNumNodePerEle(GrpEle) == 9:
            quads_o2=[]
            for quad in elelist:
                quads_o2.append([ p.id for p in quad.p])
            quadrangle_o1,big_quadrangle_o1=quadsOrdre2_to_QuadsOrdre1(quads_o2)
            tri1, X = quads_to_tris(x,y,base,X,quadrangle_o1)
        else:
            print('\nNombre de noeud par element non gere : ',mesh.getNumNodePerEle(GrpEle))
            sys.exit(1)
            
    elements = triangle + quadrangle + big_triangle_o1 + big_quadrangle_o1

    # convert all elements into triangles 
    elements_all_tris = triangle + tri1 + triangle_o1 

    # create an unstructured triangular grid instance
    triangulation = tri.Triangulation(x, y, elements_all_tris)

    # plot the contours
    plt.tricontourf(triangulation, X)

    # plot the finite element mesh
    if printMesh: plot_fem_mesh(x, y, elements)

    # show
    plt.colorbar()
    plt.axis('equal')
    plt.show()
    
    fig = plt.figure()
    ax = plt.axes(projection ='3d')
    ax.plot_trisurf(triangulation, X)

    plt.show()

