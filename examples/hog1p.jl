using NumCME 
using Catalyst
using StaticArrays: @MVector
using Sundials: CVODE_BDF

const r1 = 6.1e-3
const r2 = 6.9e-3
const η = 5.9
const Ahog = 9.3e9
const Mhog = 2.2e-2
function Hog1p(t)
    u = (1.0 - exp(-r1 * t)) * exp(-r2 * t)
    signal = Ahog * (u / (1.0 + u / Mhog))^η
    return signal 
end

@parameters begin
    k01
    k10
    a
    k12
    k21
    k23
    k32
    λ0
    λ1
    λ2
    λ3
    γnuc
    ktrans
    γcyt
end
rn = @reaction_network begin     
    k01, G0 --> G1
    max(0, k10 - a*Hog1p(t)), G1 --> G0 
    k12, G1 --> G2 
    k21, G2 --> G1 
    k23, G2 --> G3
    k32, G3 --> G2
    λ0, G0 --> G0 + RNAnuc 
    λ1, G1 --> G1 + RNAnuc
    λ2, G2 --> G2 + RNAnuc 
    λ3, G3 --> G3 + RNAnuc 
    γnuc, RNAnuc --> ∅
    ktrans, RNAnuc --> RNAcyt 
    γcyt, RNAcyt --> ∅
end 

param_values = Dict([
k01=> 2.6e-3,
k10=> 1.9e01,
a=> 0.0,
k12=> 7.63e-3,
k21=> 1.2e-2,
k23=> 4e-3,
k32=> 3.1e-3,
λ0=> 5.9e-4,
λ1=> 1.7e-1,
λ2=> 1.0,
λ3=> 3e-2,
ktrans=> 2.6e-1,
γnuc=> 2.2e-6,
γcyt=> 8.3e-3
])
# Simulate long-time behavior before MAPK signal 
model = CmeModel(rn, collect(param_values))
p0 = FspVectorSparse([@MVector [1,0,0,0,0,0]], [1.0])
fsp = AdaptiveFspSparse(
    ode_method = CVODE_BDF(linear_solver=:GMRES),
    space_adapter = RStepAdapter(10,20,true)
)
sol = solve(
    model, 
    p0,
    (0.0, 8*3600.0),
    fsp;
    saveat = [8*3600.0],
    fsptol = 1.0e-6,
    odeatol = 1.0e-14,
    odertol = 1.0e-6,
    verbose=false
);
pend = sol[end].p 
# Simulate behavior after MAPK addition 
param_values[a] = 3.2e04
model = CmeModel(rn, collect(param_values))
sol = solve(
    model, 
    pend,
    (0.0, 3600.0),
    fsp;
    saveat = [t for t in 0.0:60.0:3600.0],
    fsptol = 1.0e-4,
    odeatol = 1.0e-14,
    odertol = 1.0e-6,
    verbose=true
);

#
using Plots 

anim = @animate for i ∈ 1:length(sol)
    pfull = sol[i].p 
    pnuc = sum(pfull, [1,2,3,4,6]) |> Array 
    pcyt = sum(pfull, [1,2,3,4,5]) |> Array 
    rnaplot = plot()
    plot!(rnaplot, 0:length(pnuc)-1, pnuc, color=:blue, label="Nuclear mRNA")
    plot!(rnaplot, 0:length(pcyt)-1, pcyt, color=:red, label="Cytoplasmic mRNA")
    ylims!(rnaplot, (0.0, 1.0))
    xlims!(rnaplot, (0, 20))
    xlabel!(rnaplot, "Molecule count")
    ylabel!(rnaplot,"Probability")
    title!(rnaplot, "t = $(sol.t[i]/60.0) minute")
end
gif(anim, "hog1p.gif", fps = 10)