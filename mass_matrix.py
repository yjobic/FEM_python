# -*- coding: utf-8 -*-
"""
TD-5, bonus: the mass matrix and the numerical integration.

    M_ij = int_Omega phi_i phi_j dx      (mass matrix)
    K_ij = int_Omega grad phi_i . grad phi_j dx   (stiffness matrix)

The integrals are computed with a Gauss quadrature. For triangles,
mesh.updateBaseQuad(ordre) uses a rule which is exact for polynomials
of degree <= ordre (see FEMlib/quadrature.py).

The script prints:
  1. the eigenvalues and the condition number of M,
  2. the condition number of K (with the Dirichlet condition, as in TD-4),
  3. the smallest values of mu such that K V = mu M V (interior nodes),
  4. the error of the solution of -Delta u + c u = f, u = g on the boundary.

@author: Yann Jobic
"""

import numpy as np
import scipy.linalg
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

#local import
from FEMlib.mesh import *
from FEMlib.assembleProblem import *
from FEMlib.errors import *

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
MeshFileName = "Meshes/square_T_2_o2.msh"   # _o1 : P1, _o2 : P2
ordre = 4                                    # quadrature exact up to this degree (1 to 4)
c = 100.                                     # coefficient of the mass term

def g(x, y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)

def f(x, y):
    return (2*np.pi**2 + c)*g(x, y)

#------------------------------------------------------------------------------
# Mesh and quadrature
#------------------------------------------------------------------------------
mesh = Mesh()
mesh.GmshToMesh(MeshFileName)
mesh.updateBaseQuad(ordre)

#------------------------------------------------------------------------------
# 1. Mass matrix: eigenvalues and condition number
#------------------------------------------------------------------------------
M = coo_matrix(assembleMatrix(mesh, 0, 1).data).toarray()
eigM = np.linalg.eigvalsh(M)            # M is symmetric: real eigenvalues
tol = 1e-12*np.max(np.abs(eigM))
print("\n========== Mass matrix M ==========")
print(f"Number of dofs           : {mesh.Npts}")
print(f"Smallest eigenvalue      : {eigM.min(): .3e}")
print(f"Largest eigenvalue       : {eigM.max(): .3e}")
print(f"Zero eigenvalues         : {np.sum(np.abs(eigM) < tol)}")
print(f"Negative eigenvalues     : {np.sum(eigM < -tol)}")
print(f"Condition number of M    : {GetConditionNumber(coo_matrix(M)):.3e}")

#------------------------------------------------------------------------------
# 2. Stiffness matrix with the Dirichlet condition (as in TD-4)
#------------------------------------------------------------------------------
t = assembleMatrix(mesh, 1, 0)
B = np.zeros(mesh.Npts)
Dirichlet(mesh, t, g, B)
print(f"Condition number of K    : {GetConditionNumber(coo_matrix(t.data)):.3e}")

#------------------------------------------------------------------------------
# 3. Generalized eigenvalues K V = mu M V, on the interior nodes
#    (the boundary nodes are removed: u = 0 on the boundary)
#------------------------------------------------------------------------------
boundary = set(int(p) for gCL in mesh.CL for p in gCL)
interior = [i for i in range(mesh.Npts) if i not in boundary]
K = coo_matrix(assembleMatrix(mesh, 1, 0).data).toarray()
Ki = K[np.ix_(interior, interior)]
Mi = M[np.ix_(interior, interior)]
mu = scipy.linalg.eigvals(Ki, Mi)
mu = np.sort(np.real(mu[np.isfinite(mu)]))
print("\n========== K V = mu M V ==========")
print(f"Smallest positive mu     : {mu[mu > 0][:3]}")
print(f"Negative mu              : {mu[mu < 0][::-1][:3]}")
print(f"Infinite mu              : {len(interior) - len(mu)}")

#------------------------------------------------------------------------------
# 4. Solution of -Delta u + c u = f, u = g on the boundary
#------------------------------------------------------------------------------
t = assembleMatrix(mesh, 1, c)
B = assembleVecForceElem(mesh, f)
Dirichlet(mesh, t, g, B)
X = spsolve(coo_matrix(t.data).tocsr(), B)
SolEx = np.array([g(p.coord[0], p.coord[1]) for p in mesh.points])
print(f"\n========== Solution, c = {c} ==========")
print(f"Max error at the nodes   : {np.max(np.abs(X - SolEx)):.3e}")
