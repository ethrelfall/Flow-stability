# this is example from https://www.firedrakeproject.org/demos/rayleigh-benard.py
# with small changes by Ed Threlfall, January 2023

# solves stationary state of vertical convection problem
# attempts to do eigen-analysis of linearized perturbation d(system)/dt = - lambda*system

# transition to turbulence should happen around Ra = 3.4e5 (see paper mentioned below) so I expect to see 
# unstable growing or oscillatory perturbations around there
# ... but I don't!  What's wrong?


from firedrake import *

from firedrake.petsc import PETSc
try:
    from slepc4py import SLEPc
except ImportError:
    import sys
    warning("Unable to import SLEPc")
    sys.exit(0)

import math

# mesh as in paper "Assessment of reduced-order modelling strategies for convectie hear transfer"
# Victor Zucatti et al
M=Mesh("zucatti_paper.msh")

V = VectorFunctionSpace(M, "CG", 3)  # CG2, CG1 for pressure is Taylor-Hood
W = FunctionSpace(M, "CG", 2)
Q = FunctionSpace(M, "CG", 3)
Z = V * W * Q

upT = Function(Z)
u, p, T = split(upT)
v, q, S = TestFunctions(Z)

x, y = SpatialCoordinate(M)

Ra = Constant(1e3)
Pr = Constant(0.71)

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
    DirichletBC(Z.sub(0), Constant((0, 0)), (11, 12, 13, 14)),
    DirichletBC(Z.sub(2), Constant(1.0), (14,)),
    DirichletBC(Z.sub(2), Constant(0.0), (12,))
]

# Like Navier-Stokes, the pressure is only defined up to a constant.::

nullspace = MixedVectorSpaceBasis(
    Z, [Z.sub(0), VectorSpaceBasis(constant=True), Z.sub(2)])


# First off, we'll solve the full system using a direct solver.

from firedrake.petsc import PETSc

#continue to near Ra=3.4E5 solution then do detailed scan of eigenvalue behaviour in the vicinity

for i in range (1,12):
   Ra.assign(pow(10,(3+(i-1)/4)))
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

   print("run "+str(3+(i-1)/4))

#quit()  #TRIALCODE

f = open("zucatti_outputs.txt", "w+")
f.write("Ra, lamR, lamI\n")

for i in range (1,100):
   Ra.assign(pow(10,5.5)+i*5000)
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

   #heat flux calc - disabled
   #normL = Function(V)
   #normL = Constant((-1.0,0.0))
   #fluxL = assemble(inner(normL, grad(T))*ds(14))  
   #normR = Function(V)
   #normR = Constant((1.0,0.0))
   #fluxR = assemble(inner(normR, grad(T))*ds(12))   

   print("finished run "+str(i)+" with Ra "+str(pow(10,5.5)+i*5000))  #printing str(Ra) prints out w_3 - why?

   # do output here if desired - visualize base flow
   #u, p, T = upT.split()
   #u.rename("Velocity")
   #p.rename("Pressure")
   #T.rename("Temperature")
   #File("zucatti.pvd").write(u, p, T)

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
       + (T1*S1)*dx
   )

  #BCs zero Dirichlet on bdy where Dirichlet applies
   bcs1 = [
       DirichletBC(Z.sub(0), Constant((0, 0)), (11, 12, 13, 14)),
       DirichletBC(Z.sub(2), Constant(0.0), (14,)),
       DirichletBC(Z.sub(2), Constant(0.0), (12,))
   ]

   petsc_F1 = assemble(F1, bcs=bcs1).M.handle  #BCs here because this is wh. int by parts was applied
   petsc_m1 = assemble(m1).M.handle

   num_eigenvalues = 1

   opts = PETSc.Options()
   opts.setValue("eps_gen_non_hermitian", None)
   opts.setValue("st_pc_factor_shift_type", "NONZERO")
   opts.setValue("eps_type", "krylovschur")
   opts.setValue("eps_smallest_real", None) #real negative eigenvalues should be instability modes
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
   f.write(str(pow(10,5.5)+i*5000)+", "+str(lam1.real)+", "+str(lam1.imag)+"\n")

f.close()

quit()  #TRIALCODE

#try to output the eigenfunction (clunky, not figured out how to do properly)
import numpy as np
npa1 = np.array(vr.getSize())
npa1 = vr
cntV = V.dof_dset.size
npa2 = npa1[0:cntV]
npa3 = npa1[cntV:2*cntV]

cntW = W.dof_dset.size
npa4 = npa1[2*cntV:2*cntV+cntW]
npa5 = npa1[2*cntV+cntW:2*cntV+2*cntW]


Vcomp = FunctionSpace(M, "CG", 2)
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



