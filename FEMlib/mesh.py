# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

import sys
import numpy as np
#To install: conda install -c conda-forge gmsh python-gmsh
import gmsh
#Doc API gmsh : https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_12_1/api/gmsh.py

#local import
from FEMlib.basis import *
from FEMlib.quadrature import *


class Point:
    def __init__(self,monid, x,y):
        #coordinates
        self.coord = [x,y]
        #indice DOF, degres de liberte
        self.id=monid
    def __str__(self):
        return "id : "+str(self.id)+" coord ("+str(self.coord[0])+","+str(self.coord[1])+")"


class Ele:
    def __init__(self, dim, p, physical_tag, monid, vol):
        self.dim = dim
        assert dim==2,"On en gere pas encore le 3D, dim doit valoir 2"
        #3 points de la classe Point definissant un triangle
        self.p=p
        #surement verifier si on a bien 3 "Point" dans un array...
        self.physical_tag=physical_tag
        #ici c'est juste un indice
        self.id=monid
        #Volume de l'element
        self.vol=vol
        self.detJ=2.*vol
    def __str__(self):
        string='id : '+str(self.id)
        string+='\ntag : '+str(self.physical_tag)
        string+='\nNpts : '+str(len(self.p))
        for i in range(len(self.p)):
            string+='\npt '+str(self.p[i])
        string+='\nVol '+str(self.vol)
        return string


class Mesh:
    def __init__(self):
        self.Npts=np.int64(0)
        self.NbGroupEle=np.int32(0)
        self.infoGroups=[]
        self.points = []
        self.listElesType = []
        self.CL=[]
        #Une base par type d'elements
        self.base=[]
        #Une quadrature par type d'elements
        self.quad=[]
        
    def GmshToMesh(self, file):
        debug=0
        gmsh.initialize(sys.argv)
        gmsh.merge(file)

        #On charge les noeuds
        Nodes=gmsh.model.mesh.getNodes(dim=-1, tag=-1)
        self.Npts = len(Nodes[0])
        self.points=[ '?'  for j in range(self.Npts) ]
        print("Nombre de Noeuds : ",self.Npts)
        for i in range(self.Npts):
            self.points[i]=Point(i,Nodes[1][i*3+0],Nodes[1][i*3+1])
        
        #On charge les elements
        #getPhysicalGroups : entities are returned as a vector of (dim, tag) pairs.
        NbTotEleType=len(gmsh.model.getPhysicalGroups(dim=2))
        print("Number of physical tags of dim 2 (surfaces) : ",NbTotEleType)
        
        for i in range(NbTotEleType):
            numeroSurf=gmsh.model.getPhysicalGroups(dim=2)[i][1]
            if (debug == 1): print("\nnumSurf "+str(i)+" a pour val : ",numeroSurf)
            #On recup le num du type d'ele du groupe physique
            montag=gmsh.model.getEntitiesForPhysicalGroup(dim=2, tag=numeroSurf)[0]
            typeEle=gmsh.model.mesh.getElements(dim=2,tag=montag)[0][0]
            if (debug == 1): print('Type de l element : ',typeEle)
            elemlist=gmsh.model.mesh.getElements(dim=2,tag=montag)[1][0]
            if (debug == 1):  print('Liste des elements du groupe : ',elemlist)
            Connec=gmsh.model.mesh.getElements(dim=2,tag=montag)[2][0]
            #infos relative au type d'element
            print("Numero du type de l'element surfacique : ",typeEle,"de type "+gmsh.model.mesh.getElementProperties(typeEle)[0])
            NbNodesEle = gmsh.model.mesh.getElementProperties(typeEle)[3]
            self.infoGroups.append([NbNodesEle,numeroSurf])
            if (debug == 1): print("Nombre de noeuds par element : ",NbNodesEle)
            print("Ordre de l'element : ",gmsh.model.mesh.getElementProperties(typeEle)[2])
            
            #On reordonne la numerotation des noeuds a l'interieur de l'element
            reorder=[]
            if NbNodesEle == 3:
                reorder.append(0); reorder.append(1); reorder.append(2)
            elif NbNodesEle == 6:
                reorder.append(0); reorder.append(3); reorder.append(1);
                reorder.append(4); reorder.append(2); reorder.append(5);
            elif NbNodesEle == 4:
                reorder.append(0); reorder.append(1); reorder.append(2);
                reorder.append(3);
            elif NbNodesEle == 9:
                reorder.append(0); reorder.append(4); reorder.append(1); reorder.append(5);
                reorder.append(2); reorder.append(6); reorder.append(3); reorder.append(7);
                reorder.append(8);
            else:
                print('\n definir le mapping reorder dans GmshToMesh! \n')
                exit(1)

            listEles=[ '?'  for j in range(len(elemlist)) ]
            for j in range(len(elemlist)):
                connecEle=[ '?'  for j in range(NbNodesEle) ]
                for i in range(NbNodesEle):
                    Node=gmsh.model.mesh.getElements(dim=2,tag=montag)[2][0][NbNodesEle*j+reorder[i]]
                    #Attention : en gmsh les noeuds commencent à 1
                    N=np.int64(Node-1)
                    connecEle[i]=self.points[N]
                if (debug == 1): print(s)
                vol=gmsh.model.mesh.getElementQualities([elemlist[j]],qualityName="volume",task=0, numTasks=1)[0]
                ele=Ele(2,connecEle,numeroSurf,j,vol)
                listEles[j]=ele
            self.listElesType.append(listEles)
            
        #On fait un peu pareil pour les conditions limites
        self.numCL=len(gmsh.model.getPhysicalGroups(dim=1))
        for i in range(self.numCL):
            #On recup le tag de la condition limite
            numeroTag=gmsh.model.getPhysicalGroups(dim=1)[i][1]
            self.CL.append(gmsh.model.mesh.getNodesForPhysicalGroup(dim=1, tag=numeroTag)[0]-1)

        # Finalize GMSH
        print("")
        gmsh.finalize()
        self.NbGroupEle=len(self.infoGroups)
            
    def updateBaseQuad(self,ordre):
        for gp in range(self.NbGroupEle):
            NumNodePerEle=self.getNumNodePerEle(gp)
            #switch case only working for python>=3.10
            if NumNodePerEle == 3:
                self.base.append(LagrangeT1(2,NumNodePerEle))
                self.quad.append(QuadratureT(2,ordre))
                print('On ajoute au groupe ',gp,'la base ',self.base[gp].name,'\n')
            elif NumNodePerEle == 6:
                self.base.append(LagrangeT2(2,NumNodePerEle))
                self.quad.append(QuadratureT(2,ordre))
                print('On ajoute au groupe ',gp,'la base ',self.base[gp].name,'\n')
            elif NumNodePerEle == 4:
                self.base.append(LagrangeQ1(2,NumNodePerEle))
                self.quad.append(QuadratureQ(2,ordre))
                print('On ajoute au groupe ',gp,'la base ',self.base[gp].name,'\n')
            elif self.getNumNodePerEle(gp) == 9:
                self.base.append(LagrangeQ2(2,NumNodePerEle))
                self.quad.append(QuadratureQ(2,ordre))
            else:
                print('\nNombre de Noeud par element non connu : ',self.getNumNodePerEle(gp))
                sys.exit(1)

    def printNodes(self):
        print("Nombre de Noeuds : ",self.Npts)
        for i in range(self.Npts):
            print(self.points[i])
    
    def printConnec(self):
        print("Nombre de groupe d'elements : ",self.NbGroupEle)
        for i,gp in zip(range(self.NbGroupEle),self.infoGroups):
            nNodeParEle=gp[0]
            numGroup=gp[1]
            print("Group ",str(i)," NumGroup ",str(numGroup)," Type ele ",nNodeParEle)
            elemlist=self.listElesType[i]
            Nele=len(elemlist)
            #NbPtParEle=len(Connec[0])
            for j in range(Nele):
                s="Ele "+str(j)+" a pour connecte les pts : "
                for pt in elemlist[j].p:
                    s+=" "+str(pt.id)
                print(s)
           
    def printInfoMesh(self):
        print("Nombre de groupe d'elements : ",self.NbGroupEle)
        for i,gp in zip(range(self.NbGroupEle),self.infoGroups):
            print("group "+str(i)+", tag du groupe : "+str(gp[1])+", nombre de noeuds par ele : ",str(gp[0]),'qui a ',str(len(self.listElesType[i])),' elements')
            #print('groupe num ',i,'Num gp[0],gp[1])

    def getNumNodePerEle(self, groupId):
        return self.infoGroups[groupId][0]
    
    def getTagGroupe(self,groupId):
        return self.infoGroups[groupId][1]

    def loc2Glob(self, ele, i):
        return ele.p[i].id
    
    def getTotEle(self):
        tot=0
        for gp in range(self.NbGroupEle):
            tot+=len(self.listElesType[gp])
        return tot
        
