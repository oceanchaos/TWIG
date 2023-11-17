module m

import Models: Model
using ParametricModels
using LinearAlgebra
using SpecialFunctions
import Sundials

@parameterspace mutable struct Bif
   y0		= .5, identity
   r		= 0, identity
   alpha1	= 0, identity
   alpha2	= 0, identity
   alpha3	= 0, identity
   alpha4	= 0, identity
   alpha5	= 0, identity
end

### Initial condition function
function ic(ps::Bif{T}) where T<:Real
   return T[ps.y0]
end

### define right hand side function
function rhs(ps::Bif{T}, t, y, dy) where T<:Real
	dy[1]=ps.r*y[1]-y[1]^3+ps.alpha1*y[1]^4+ps.alpha2*y[1]^5+ps.alpha3*y[1]^6+ps.alpha4*y[1]^7+ps.alpha5*y[1]^8	
    nothing
end

### define observation function
function obs(ps::Bif{T}, t, y) where T<:Real
   return y
end

jacobians = Matrix[]
coords = Any[]
lambdas = Any[]
eig=Any[]
lsing=Any[]

tmaxs = 10.0 .^ (-2:.1:3)
for (i, tmax)  = enumerate(tmaxs)
	@info "i=$i, tmax=$(tmax)"
	t=range(0, stop = tmax, length = 101) |>collect
	f(ps::Bif{T}, t) where T <: Real = solve_ode(ps, ic, rhs, obs, t[2:end], Sundials.CVODE_BDF(); abstol = 1e-14, reltol =1e-14)

	pmodel = PModel(Bif, parameter_transforms, f, OLSData("bifdata", ModelArgs(t), f(Bif(), t)))
	model = Model(pmodel)
	x = xvalues(pmodel)

	j = model.jacobian(x)
	push!(jacobians, j)
	push!(coords, pmodel(x)[1])
	u,s,v=svd(j)
	push!(lambdas, s)
	push!(eig,v)
	push!(lsing,u)
 end
end

using DelimitedFiles 
eigs=hcat(m.eig...)
lam=hcat(m.lambdas...)
writedlm("eigvec_pitchf.txt",eigs)
writedlm("eigval_pitchf.txt",lam)
writedlm("tmaxs.txt",m.tmaxs)
writedlm("lsing_pitchf.txt",hcat(m.lsing...))