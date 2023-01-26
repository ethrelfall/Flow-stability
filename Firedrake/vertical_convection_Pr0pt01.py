# this is example from https://www.firedrakeproject.org/demos/rayleigh-benard.py
# with small changes by Ed Threlfall, January 2023


from firedrake import *

from firedrake.petsc import PETSc
try:
    from slepc4py import SLEPc
except ImportError:
    import sys
    warning("Unable to import SLEPc")
    sys.exit(0)

import math

M = Mesh("vertical_convection_40.msh")  # 40*40, stretched in x and y warp factor w=12
#M = Mesh("vertical_convection_80.msh") #80*80, warp factor 12 

V = VectorFunctionSpace(M, "CG", 6)
W = FunctionSpace(M, "CG", 5)
Q = FunctionSpace(M, "CG", 6)
Z = V * W * Q

upT = Function(Z)
u, p, T = split(upT)
v, q, S = TestFunctions(Z)

x, y = SpatialCoordinate(M)

Ra = Constant(1e5)
Pr = Constant(0.01)  # this is like a liquid metal

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

f = open("convection_outputsPr0pt01.txt", "w+")
f.write("Ra, Nu, lamR, lamI\n")

for i in range (1,48):
   Ra.assign(pow(10,(5+(i-1)/40)))
   print("starting run "+str(i)+" with Ra "+str((5+(i-1)/40)))
   try:
      solve(F == 0, upT, bcs=bcs, nullspace=nullspace,
            solver_parameters={"mat_type": "aij",
                               "snes_monitor": None,
                               "ksp_type": "gmres",
                               "pc_type": "lu",
                               "pc_factor_mat_solver_type": "mumps"})

   except:
      print("something bad happened on run "+str(i))

   normL = Function(V)
   normL = Constant((-1.0,0.0))
   fluxL = assemble(inner(normL, grad(T))*ds(14))  
   normR = Function(V)
   normR = Constant((1.0,0.0))
   fluxR = assemble(inner(normR, grad(T))*ds(12))   
   print("finished run "+str(i)+" with Ra "+str(Ra))

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

  #BCs zero on Dirichlet bdy
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
   opts.setValue("eps_largest_imaginary", None)
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
   lam1 = es.getEigenpair(0, vr, vi)  #choose which mode
   print("eigenvalue 1 value:")
   print(lam1)

   f.write(str(5+(i-1)/40)+", "+str(math.log10(abs(0.5*(fluxL-fluxR))))+", "+str(lam1.real)+", "+str(lam1.imag)+"\n")

   #TRIALCODE
   #u, p, T = upT.split()
   #u.rename("Velocity")
   #p.rename("Pressure")
   #T.rename("Temperature")
   #File("benard_test.pvd").write(u,p,T)

f.close()

#do output here if desired
#u, p, T = upT.split()
#u.rename("Velocity")
#p.rename("Pressure")
#T.rename("Temperature")
#File("benard_test.pvd").write(u, p, T)

quit()

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


Vcomp = FunctionSpace(M, "CG", 4)
eigenmode = Function(Vcomp)
eigenmode.vector()[:] = npa2
eigenmode2 = Function(Vcomp)
eigenmode2.vector()[:] = npa3
File("benard_eigenmode_u.pvd").write(eigenmode, eigenmode2)

eigenmode3 = Function(W)
eigenmode3.vector()[:] = npa4
File("benard_eigenmode_p.pvd").write(eigenmode3)

eigenmode4 = Function(Q)
eigenmode4.vector()[:] = npa5
File("benard_eigenmode_T.pvd").write(eigenmode4)




print('finished - quitting.')
quit()


