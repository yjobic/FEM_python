import ufl
from dolfinx import mesh,fem,default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import numpy as np

def g(x,y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)

def f(x,y):
    return 2*np.pi**2*g(x,y)

domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8, mesh.CellType.quadrilateral)
V = fem.functionspace(domain, ("Lagrange", 1))

uD = fem.Function(V)
uD.interpolate(lambda x: g(x[0],x[1]) )

tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)

boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
bc = fem.dirichletbc(uD, boundary_dofs)


#Trial and test function
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

#Variational problem
a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = g * v * ufl.dx

problem = LinearProblem(
    a,
    L,
    bcs=[bc],
    petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    petsc_options_prefix="Poisson",
)
uh = problem.solve()

V2 = fem.functionspace(domain, ("Lagrange", 2))
ue = fem.Function(V2, name="u_exact")
ue.interpolate(lambda x: g(x[0],x[1]) )

L2_error = fem.form(ufl.inner(uh - ue, uh - ue) * ufl.dx)
error_local = fem.assemble_scalar(L2_error)
error_L2 = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))

local_L_infty = np.max(np.abs(uD.x.array - uh.x.array))
L_infty = domain.comm.allreduce(local_L_infty, op=MPI.MAX)
if domain.comm.rank == 0:  # Only print the error on one process
    print(f"Error_L2 : {error_L2:.2e}")
    print(f"Error_max : {L_infty:.2e}")


