# -*- coding: utf-8 -*-
"""
Linear elasticity with FEniCSx: 2D cantilever beam under its own weight
(plane stress). The beam is clamped on the left side (x=0).

    -div(sigma(u)) = f    in Omega = ]0,L[ x ]0,H[
                 u = 0    on x = 0
     sigma(u) . n  = 0    on the other sides

@author: Yann Jobic
"""

import ufl
from dolfinx import mesh, fem, plot, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import numpy as np
import pyvista

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
Lx, H = 1.0, 0.1            # length and height of the beam (m)
E, nu = 210e9, 0.3          # Young modulus (Pa), Poisson ratio (steel)
rho, grav = 7800., 9.81     # density (kg/m^3), gravity (m/s^2)
degree = 1                  # degree of the Lagrange polynomials (1 or 2)
ny = 4                      # number of cells in the thickness
nx = 10*ny                  # cells of aspect ratio 1

# Lame coefficients, plane stress
mu = E/(2*(1 + nu))
lmbda = E*nu/(1 - nu**2)

#------------------------------------------------------------------------------
# Mesh and vector function space (2 components per node)
#------------------------------------------------------------------------------
domain = mesh.create_rectangle(MPI.COMM_WORLD, [[0., 0.], [Lx, H]], [nx, ny],
                               mesh.CellType.triangle)
gdim = domain.geometry.dim
V = fem.functionspace(domain, ("Lagrange", degree, (gdim,)))

#------------------------------------------------------------------------------
# Clamped boundary: u = (0,0) on x = 0
#------------------------------------------------------------------------------
clamped_dofs = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x[0], 0.))
bc = fem.dirichletbc(np.zeros(gdim, dtype=default_scalar_type), clamped_dofs, V)

#------------------------------------------------------------------------------
# Weak formulation
#------------------------------------------------------------------------------
def epsilon(u):
    return ufl.sym(ufl.grad(u))

def sigma(u):
    return lmbda*ufl.tr(epsilon(u))*ufl.Identity(gdim) + 2*mu*epsilon(u)

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, default_scalar_type((0., -rho*grav)))   # weight

#==============================================================================
# TO DO: write the weak formulation
#     a(u,v) = int_Omega sigma(u) : epsilon(v) dx
#     L(v)   = int_Omega f . v dx
# Useful ufl functions: ufl.inner (A : B), ufl.dot (a . b), ufl.dx
# Compare with heat.py!
#==============================================================================
a = None
L = None

if a is None or L is None:
    raise NotImplementedError("Write the weak formulation (a and L) first!")

problem = LinearProblem(a, L, bcs=[bc],
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
                        petsc_options_prefix="Elasticity")
uh = problem.solve()
uh.name = "displacement"

#------------------------------------------------------------------------------
# Tip deflection vs beam theories
#------------------------------------------------------------------------------
q = rho*grav*H                  # load per unit length (unit thickness)
I = H**3/12                     # second moment of area (unit thickness)
delta_EB = q*Lx**4/(8*E*I)      # Euler-Bernoulli
kappa = 5./6                    # shear correction factor
delta_T = delta_EB + q*Lx**2/(2*kappa*mu*H)   # Timoshenko

# Vertical displacement of the point (Lx, H/2)
coords = V.tabulate_dof_coordinates()
tip = np.argmin(np.linalg.norm(coords[:, :2] - [Lx, H/2], axis=1))
delta_h = -uh.x.array.reshape(-1, gdim)[tip, 1]

num_dofs = V.dofmap.index_map.size_global*V.dofmap.index_map_bs
print(f"P{degree}, {nx}x{ny} cells, {num_dofs} dofs")
print(f"delta FEM : {delta_h:.4e} m")
print(f"delta Euler-Bernoulli : {delta_EB:.4e} m  (ratio {delta_h/delta_EB:.3f})")
print(f"delta Timoshenko : {delta_T:.4e} m  (ratio {delta_h/delta_T:.3f})")

#------------------------------------------------------------------------------
# Visualization with pyvista: deformed shape (amplified)
#------------------------------------------------------------------------------
topology, cell_types, geometry = plot.vtk_mesh(V)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
values = np.zeros((geometry.shape[0], 3))       # pyvista needs 3D vectors
values[:, :gdim] = uh.x.array.real.reshape(-1, gdim)
grid.point_data["u"] = values
grid.set_active_vectors("u")
warped = grid.warp_by_vector("u", factor=0.2*Lx/delta_EB)   # deflection shown as 20% of L

plotter = pyvista.Plotter(window_size=(1200, 400))
plotter.add_mesh(grid, style="wireframe", color="gray")
plotter.add_mesh(warped, scalars="u", show_edges=True,
                 scalar_bar_args={"title": "|u| (m)", "fmt": "%.1e", "n_labels": 3})
plotter.view_xy()
plotter.camera.zoom(1.5)
if pyvista.OFF_SCREEN:
    plotter.screenshot(f"beam_P{degree}.png")
else:
    plotter.show()
