using Random
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Visualization
using GridapGmsh
using Printf

path = "/mnt/c/Users/Alvaro Baillet Boliv/vtu_files/"

# For testing only
abstract type Estimator end

struct RandomEst <: Estimator
  function RandomEst(seed)
    Random.seed!(seed)
    new()
  end
end

function choose_comparison(
    model::DiscreteModel,
    hx::Float64,
    hy::Float64
)
    #do stuff
    coords = get_cell_coordinates(model)

    for c in coords
        mid = (c[1] + c[2] + c[3])/3  
        mid_adj = mid - -1 # assuming low left point is -1, -1
        quad = (Integer(mid_adj[1] ÷ hx), Integer(mid_adj[2] ÷ hy)) # quad holds which rectangle contains element

        # Initial mesh has downwards division ie: |\|. So to check which triangle:
        x1 = quad[1] * hx + -1  # the + 0 is temporary, assuming the domain starts at x = -1
        y1 = (quad[2] + 1) * hy + -1 # + 0 assumes domain starts at y = -1
        x2 = (quad[1] + 1) * hx + -1
        y2 = quad[2] * hy + -1
        slope = (y2 - y1)/(x2 - x1)
        
        if (mid[2] < slope * (mid[1] - x1) + y1)
            idx = (quad[1], quad[2], 0) # idx contains which triangle contains element
        else
            idx = (quad[1], quad[2], 1) # idx contains which triangle contains element
        end
        print(idx,"\t")
        print(c[1]," ", c[2], " ", c[3],"\n")
    end
end

function make_nvb_levels(
  model::DiscreteModel,
  Nsteps::Integer,
  θ::AbstractFloat,
  est::Estimator,
)
  model_refs = Vector{DiscreteModel}(undef, Nsteps)
  cell_map = get_cell_map(get_triangulation(model))
  ncells = length(cell_map)
  print("1: "*string(ncells)*"\n")
  η_arr = compute_estimator(est, ncells)
  model_refs[1], buffer = newest_vertex_bisection(model, η_arr; θ = θ)
  buffer = deepcopy(buffer)
  for i = 1:(Nsteps - 1)
    cell_map = get_cell_map(get_triangulation(model_refs[i]))
    ncells = length(cell_map)
    print(string(i+1)*": "*string(ncells)*"\n")
    η_arr = compute_estimator(est, ncells)
    model_refs[i + 1], buffer =
      newest_vertex_bisection(model_refs[i], buffer, η_arr; θ = θ)
  # Necessary to have each level stored seperately
  buffer = deepcopy(buffer)
  end
  model_refs
end

compute_estimator(est::RandomEst, ncells) = rand(ncells)


# Nonuniform refinement.
domain = (0, 1, 0, 1)
partition = (2, 2) # Initial partition
model = simplexify(CartesianDiscreteModel(domain, partition))
writevtk(model, path*"initial_model") 

model_coord_lowleft = get_cell_coordinates(model)[1]

hx = abs((model_coord_lowleft[2] - model_coord_lowleft[1])[1])
hy = abs((model_coord_lowleft[3] - model_coord_lowleft[1])[2])

Nsteps = 10
seed = 5
est = RandomEst(seed)
θ = 0.5
nonuniform_write_to_vtk = true


@time model_refs = make_nvb_levels(model, Nsteps, θ, est)
if nonuniform_write_to_vtk
 for (n, model_ref) in enumerate(model_refs)
   trian_ref = get_triangulation(model_ref)
   writevtk(trian_ref, path*"nonuniform$(string(n, pad=2))")
 end
end

test_model = model_refs[1]


choose_comparison(test_model, hx, hy) 