# -*- coding: utf-8 -*-
"""
Solve with FEniCSx the same problem as with FEM_python:

    -Delta u = f   in Omega = ]0,1[^2
           u = g   on the boundary

with the manufactured solution g(x,y) = sin(pi x) sin(pi y),
hence f = (2 pi^2) g.
c = 0 gives the Poisson (heat) equation of TD-4.

@author: Yann Jobic
"""

import ufl
from dolfinx import mesh, fem, plot, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
import numpy as np
import pyvista

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
c = 0.              # coefficient of the mass term (c=0 : Poisson)
degree = 1          # degree of the Lagrange polynomials (1 or 2)
MeshFileName = None # e.g. "../Meshes/square_T_3_o1.msh" to use the Gmsh meshes of FEM_python
N = 8               # number of cells in each direction, if MeshFileName is None

#------------------------------------------------------------------------------
# Exact solution and source term
#------------------------------------------------------------------------------
# The module "m" is a parameter, so that the same function can be evaluated:
#  - with numpy (m=np) on arrays of points, to interpolate a fem.Function,
#  - with ufl (m=ufl) on the symbolic coordinates, to write the variational form.
def g(x, y, m=np):
    return m.sin(m.pi*x)*m.sin(m.pi*y)

def f(x, y, m=np):
    return (2*m.pi**2)*g(x, y, m)

#------------------------------------------------------------------------------
# Mesh
#------------------------------------------------------------------------------
if MeshFileName is None:
    # Structured mesh of the unit square (quadrilaterals or triangles)
    domain = mesh.create_unit_square(MPI.COMM_WORLD, N, N, mesh.CellType.quadrilateral)
else:
    # Same Gmsh mesh as in FEM_python
    domain = gmshio.read_from_msh(MeshFileName, MPI.COMM_WORLD, 0, gdim=2).mesh

# Finite element space: continuous Lagrange polynomials of degree "degree"
V = fem.functionspace(domain, ("Lagrange", degree))

#------------------------------------------------------------------------------
# Dirichlet boundary condition
#------------------------------------------------------------------------------
# uD is the finite element function interpolating g (nodal values)
uD = fem.Function(V)
uD.interpolate(lambda x: g(x[0], x[1]))

# Find the boundary facets (edges in 2D) and the dofs located on them
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)
boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
bc = fem.dirichletbc(uD, boundary_dofs)

#------------------------------------------------------------------------------
# Variational problem
#------------------------------------------------------------------------------
# Trial and test functions
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# Symbolic coordinates (x[0], x[1]) of the domain: f is written as a ufl
# expression, evaluated exactly at the quadrature points during assembly.
x = ufl.SpatialCoordinate(domain)
f_ufl = f(x[0], x[1], ufl)

# Bilinear form: stiffness term + c * mass term
c_ = fem.Constant(domain, default_scalar_type(c))
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
# Linear form
L = f_ufl * v * ufl.dx

# Assembly and resolution with a direct solver (LU)
problem = LinearProblem(
    a,
    L,
    bcs=[bc],
    petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    petsc_options_prefix="Poisson",
)
uh = problem.solve()
uh.name = "u_h"

#------------------------------------------------------------------------------
# Errors
#------------------------------------------------------------------------------
# L2 error: the exact solution is also written as a ufl expression, so the
# error is not polluted by an interpolation of g.
# The quadrature degree is increased because g is not a polynomial.
g_ufl = g(x[0], x[1], ufl)
L2_error = fem.form((uh - g_ufl)**2 * ufl.dx(metadata={"quadrature_degree": 2*degree + 3}))
error_local = fem.assemble_scalar(L2_error)
error_L2 = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))

# Maximum error at the nodes (uD contains the exact nodal values)
local_L_infty = np.max(np.abs(uD.x.array - uh.x.array))
L_infty = domain.comm.allreduce(local_L_infty, op=MPI.MAX)

num_cells = domain.topology.index_map(tdim).size_global
num_dofs = V.dofmap.index_map.size_global
if domain.comm.rank == 0:  # Only print the error on one process
    print(f"Number of cells : {num_cells}")
    print(f"Number of dofs : {num_dofs}")
    print(f"Error_L2 : {error_L2:.2e}")
    print(f"Error_max : {L_infty:.2e}")

#------------------------------------------------------------------------------
# Visualization with pyvista
#------------------------------------------------------------------------------
# Convert the function space into a VTK unstructured grid: one point per dof
topology, cell_types, geometry = plot.vtk_mesh(V)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
grid.point_data["u_h"] = uh.x.array.real
grid.point_data["error"] = np.abs(uh.x.array.real - uD.x.array.real)
grid.set_active_scalars("u_h")

plotter = pyvista.Plotter(shape=(1, 3), window_size=(1500, 500))

# Left: solution on the mesh
plotter.subplot(0, 0)
plotter.add_text("u_h", font_size=10)
plotter.add_mesh(grid.copy(), scalars="u_h", show_edges=True,
                 scalar_bar_args={"title": "u_h", "fmt": "%.2f"})
plotter.view_xy()

# Middle: solution as a surface z = u_h(x,y)
plotter.subplot(0, 1)
plotter.add_text("u_h (warped)", font_size=10)
warped = grid.warp_by_scalar("u_h", factor=0.5)
plotter.add_mesh(warped, scalars="u_h", show_scalar_bar=False)

# Right: nodal error
plotter.subplot(0, 2)
plotter.add_text("|u_h - g| at the nodes", font_size=10)
plotter.add_mesh(grid.copy(), scalars="error", show_edges=True, cmap="plasma",
                 scalar_bar_args={"title": "error", "fmt": "%.1e", "n_labels": 3})
plotter.view_xy()

if pyvista.OFF_SCREEN:
    plotter.screenshot("solution.png")
else:
    plotter.show()
