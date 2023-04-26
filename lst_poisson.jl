using Gridap
using GridapGmsh

path = "/mnt/c/Users/Alvaro Baillet Boliv/vtu_files/" # this is only for Alvaro, to move the file to windows
# to run the code delete the path from lines 15 and 46 (where the writevtk func is)



# Geometry
domain = (0,1,0,1) # lines 10-12 define the mesh
partition = (4,4)
model = simplexify(CartesianDiscreteModel(domain,partition))
print("Done model\n")
Ω = Triangulation(model)
writevtk(model, path*"model")

print("Done geometry \n")

# Interpolation
reffe = ReferenceFE(lagrangian, Float64, 1) #order 1 Finite Elements 
V = FESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
g(x) = 2.0 # dirichlet b.c.
U = TrialFESpace(V, g)

print("Done Interpolation \n")

# Quadrature
dΩ = Measure(Ω, 2) #quadrature order 2

# Weak form
f(x) = 6*x[1] + 4*x[2] # source function
a(u, v) = ∫(∇(u)⋅∇(v)) * dΩ
l(v) = ∫(v* f) *dΩ

print("Done weak \n")

# assembly
op = AffineFEOperator(a, l, U, V)
print("Done Assembly\n")
solver = FESolver(LUSolver())

uh = solve(solver, op)

print("Done Solve \n")

writevtk(Ω, path*"results_own2", cellfields=["uh"=>uh])



