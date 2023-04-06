### Read the Compass slice from SLIM website and take out the 2D slice
## Please refer to the details in the MIT license of this repository and in the license of the Compass model
## author: Ziyi Yin (ziyi.yin@gatech.edu)

using DrWatson
@quickactivate "proxy-example"
using Pkg; Pkg.instantiate()
using PyPlot
using SegyIO
using Statistics
using Polynomials
using JutulDarcyRules

include(srcdir("utils.jl"))
if ~isdir("slim.gatech.edu/data//synth/Compass")
    run(`wget -r ftp://slim.gatech.edu/data//synth/Compass`) # this might take a while
end

# get velocity
block = segy_read("slim.gatech.edu/data/synth/Compass/final_velocity_model_ieee_6m.sgy")

# original compass model is in 25m*25m*6m
n = (1911,2730,341)
d = (25f0,25f0,6f0)

sx = get_header(block, "SourceX")
sy = get_header(block, "SourceY")

#= in case want to turn it to 3D cube
v = zeros(Float32,n)
for i = 1:n[1]
    x = d[1].*(i-1)
    inds = findall(sx.==x)
    slice = block.data[:,inds[sortperm(sy[inds])]]

    v[i,:,:] = transpose(slice)/1f3
end
=#

# take a slice for now
x = 0.
inds = findall(sx.==x)
slice = block.data[:,inds[sortperm(sy[inds])]]'/1f3

factor = (3,2)                                                              # downsampling
vsmall = Float64.(1f0./downsample(1f0./slice[1:768, 182:end], (3,2)))        # 2D slice
h = 181.0 * d[end]                                                            # starting depth of the simulation model

Kh = VtoK.(vsmall)          # horizontal permeability
K = Kh * md                 # millidary
ϕ = Ktoϕ.(Kh)                # porosity

kvoverkh = 0.1              # kv/kh ratio
n = (size(K,1), 1, size(K,2))                   # model size for flow simulation
d = (Float64(d[1] * factor[1]), Float64(d[1] * factor[1] * n[1]/5), Float64(d[end] * factor[end]))  # discretization

model = jutulModel(n, d, vec(padϕ(ϕ)), K1to3(K; kvoverkh=kvoverkh), h)

## simulation time steppings
tstep = 365.25 * 12 * ones(5)       # 60 years
tot_time = sum(tstep)

## injection well location
pore_volumes = sum(ϕ[2:end-1,1:end-1] .* (vsmall[2:end-1,1:end-1].>3.5)) * prod(d)      ## reservoir pore volumes (v>3.5 is reservoir)
irate = 0.2 * pore_volumes / tot_time / 24 / 60 / 60                                    ## set injection rate for 20% storage capacity
q = jutulVWell(irate, (128 * d[1], 1 * d[2]); startz = (n[end]-18) * d[end], endz = (n[end]-16) * d[end])       ## well with rates and locations

## set up modeling operator
S = jutulModeling(model, tstep)

## simulation
mesh_ = CartesianMesh(model)            ## cartesian mesh
T(x) = log.(KtoTrans(mesh_, K1to3(exp.(x); kvoverkh=kvoverkh)))         ## convert log permeability to log transmissibility (which jutul prefers)

## log permeability
logK = log.(K)

## simulation
@time state = S(T(logK), vec(padϕ(ϕ)), q)

### plotting
extent = (0f0, (n[1]-1)*d[1], (n[end]-1)*d[end]+h, 0f0+h)
figure(figsize=(10,6))
subplot(3,1,1);
scatter(range(q.loc[1][1], stop=q.loc[1][1], length=3), h .+ range(q.startz, stop=q.endz, length=3), label="injection well", marker="x", color="black", s=30);
legend()
imshow(exp.(logK)'./md, norm=matplotlib.colors.LogNorm(vmin=100, vmax=maximum(exp.(logK)./md)), extent=extent, aspect="auto"); xlabel("X [m]"); ylabel("Z [m]"); colorbar(); title("true permeability");
subplot(3,1,2);
scatter(range(q.loc[1][1], stop=q.loc[1][1], length=3), h .+ range(q.startz, stop=q.endz, length=3), label="injection well", marker="x", color="black", s=30);
legend()
imshow(reshape(Saturations(state.states[end]), n[1], n[end])', extent=extent, aspect="auto", cmap="gnuplot"); xlabel("X [m]"); ylabel("Z [m]"); colorbar();title("saturation after $(tot_time/365.25) years");
subplot(3,1,3);
scatter(range(q.loc[1][1], stop=q.loc[1][1], length=3), h .+ range(q.startz, stop=q.endz, length=3), label="injection well", marker="x", color="black", s=30);
legend()
imshow(reshape(Pressure(state.states[end]), n[1], n[end])', extent=extent, aspect="auto"); xlabel("X [m]"); ylabel("Z [m]"); colorbar();title("pressure after $(tot_time/365.25) years");
tight_layout()
mkpath(plotsdir())
savefig(plotsdir("sat_p.png"), bbox_inches="tight", dpi=300)
