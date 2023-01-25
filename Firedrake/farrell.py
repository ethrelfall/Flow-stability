# this is example from https://www.firedrakeproject.org/demos/rayleigh-benard.py
# with small changes by Ed Threlfall, January 2023

# finds bifurcation instability in Rayleigh-Benard convection as in
# Section IV.A of
# "Bifucation analysis of two-dimensional Rayleigh-Benard convection using deflation"
# by Boulle, Dallas, and Farrell.

from firedrake import *

from firedrake.petsc import PETSc
try:
    from slepc4py import SLEPc
except ImportError:
    import sys
    warning("Unable to import SLEPc")
    sys.exit(0)

import math

N=50
M = UnitSquareMesh(50, 50)

V = VectorFunctionSpace(M, "CG", 2)  # CG2, CG1 for pressure is Taylor-Hood
W = FunctionSpace(M, "CG", 1)
Q = FunctionSpace(M, "CG", 2)
Z = V * W * Q

upT = Function(Z)
u, p, T = split(upT)
v, q, S = TestFunctions(Z)

x, y = SpatialCoordinate(M)

Ra = Constant(1e4)
Pr = Constant(1.0)

g = Constant((0, 1))  # Ed: changed sign here because the buoyancy force is UP

F = (
    Pr * inner(grad(u), grad(v))*dx
    + inner(dot(grad(u), u), v)*dx
    - inner(p, div(v))*dx
    - (Ra*Pr)*inner(T*g, v)*dx
    + inner(div(u), q)*dx
    + inner(dot(grad(T), u), S)*dx
    + inner(grad(T), grad(S))*dx
)

bcs = [
    DirichletBC(Z.sub(0), Constant((0, 0)), (1, 2, 3, 4)),
    DirichletBC(Z.sub(2), Constant(1.0), (3,)),
    DirichletBC(Z.sub(2), Constant(0.0), (4,))
]

# Like Navier-Stokes, the pressure is only defined up to a constant.::

nullspace = MixedVectorSpaceBasis(
    Z, [Z.sub(0), VectorSpaceBasis(constant=True), Z.sub(2)])


# First off, we'll solve the full system using a direct solver.

from firedrake.petsc import PETSc

f = open("farrell_outputs.txt", "w+")
f.write("Ra, lamR, lamI\n")

for i in range (1,30):
   Ra.assign(i*100)
   print("starting run.")
   try:
      solve(F == 0, upT, bcs=bcs, nullspace=nullspace,
            solver_parameters={"mat_type": "aij",
                               "snes_monitor": None,
                               "ksp_type": "gmres",
                               "pc_type": "lu",
                               "pc_factor_mat_solver_type": "mumps"})

   except:
      print("something bad happened on run "+str(i))

   print("finished run "+str(i)+" with Ra "+str(i*100))

   # do output here if desired
   u, p, T = upT.split()
   u.rename("Velocity")
   p.rename("Pressure")
   T.rename("Temperature")
   File("farrell.pvd").write(u, p, T)

   #here is the stability eigen-analysis

   upT1= TrialFunction(Z)  #must be TrialFunction not just Function for SLEPc to work!
   u1, p1, T1 = split(upT1)
   v1, q1, S1 = TestFunctions(Z)

   F1 = (
       Pr * inner(grad(u1), grad(v1))*dx
       + inner(dot(grad(u), u1), v1)*dx
       + inner(dot(grad(u1), u), v1)*dx
       - inner(p1, div(v1))*dx
       - (Ra*Pr)*inner(T1*g, v1)*dx
       + inner(div(u1), q1)*dx
       + inner(dot(grad(T), u1), S1)*dx
       + inner(dot(grad(T1), u), S1)*dx
       + inner(grad(T1), grad(S1))*dx
   )

   m1 = (
       inner(u1,v1)*dx
       +(T1*S1)*dx
   )

  #BCs zero on Dirichlet bdy
   bcs1 = [
       DirichletBC(Z.sub(0), Constant((0, 0)), (1, 2, 3, 4)),
       DirichletBC(Z.sub(2), Constant(0.0), (3,)),
       DirichletBC(Z.sub(2), Constant(0.0), (4,))
   ]

   petsc_F1 = assemble(F1, bcs=bcs1).M.handle  #BCs here because this is wh. int by parts was applied
   petsc_m1 = assemble(m1).M.handle

   num_eigenvalues = 1

   opts = PETSc.Options()
   opts.setValue("eps_gen_non_hermitian", None)
   opts.setValue("st_pc_factor_shift_type", "NONZERO")
   opts.setValue("eps_type", "krylovschur")
   opts.setValue("eps_smallest_real", None)
   opts.setValue("eps_tol", 1e-10)

   es = SLEPc.EPS().create(comm=COMM_WORLD)
   es.setDimensions(num_eigenvalues)
   es.setOperators(petsc_F1, petsc_m1)
   es.setFromOptions()

   st=es.getST()
   st.setType(SLEPc.ST.Type.SINVERT)
   es.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)

   es.solve()

   nconv = es.getConverged()
   print("number of converged eigenvalues:")
   print(nconv)   
   vr, vi = petsc_F1.getVecs()

   for i in range (0, nconv):
      lam = es.getEigenpair(i, vr, vi)  #choose which mode
      print("eigenvalue i value:")
      print(lam)


   lam1 = es.getEigenpair(0, vr, vi)
   f.write(str(i*100)+", "+str(lam1.real)+", "+str(lam1.imag)+"\n")

   #TRIALCODE
   #u, p, T = upT.split()
   #u.rename("Velocity")
   #p.rename("Pressure")
   #T.rename("Temperature")
   #File("benard_test.pvd").write(u,p,T)

f.close()

quit()

#do output here if desired
#u, p, T = upT.split()
#u.rename("Velocity")
#p.rename("Pressure")
#T.rename("Temperature")
#File("benard_test.pvd").write(u, p, T)

#try to output the eigenfunction (clunky)
import numpy as np
npa1 = np.array(vr.getSize())
npa1 = vr
cntV = V.dof_dset.size
npa2 = npa1[0:cntV]
npa3 = npa1[cntV:2*cntV]

cntW = W.dof_dset.size
npa4 = npa1[2*cntV:2*cntV+cntW]
npa5 = npa1[2*cntV+cntW:2*cntV+2*cntW]


Vcomp = FunctionSpace(M, "CG", 8)  #gotta match V ...
eigenmode = Function(Vcomp)
eigenmode.vector()[:] = npa2
eigenmode2 = Function(Vcomp)
eigenmode2.vector()[:] = npa3
File("zucatti_eigenmode_u.pvd").write(eigenmode, eigenmode2)

eigenmode3 = Function(W)
eigenmode3.vector()[:] = npa4
File("zucatti_eigenmode_p.pvd").write(eigenmode3)

eigenmode4 = Function(Q)
eigenmode4.vector()[:] = npa5
File("zucatti_eigenmode_T.pvd").write(eigenmode4)




print('finished - quitting.')
quit()


