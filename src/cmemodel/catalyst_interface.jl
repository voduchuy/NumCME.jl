function CmeModel(rn::Catalyst.ReactionSystem, parameter_values)
    species_count = Catalyst.numspecies(rn)
    reaction_count = Catalyst.numreactions(rn)
    parameter_count = Catalyst.numreactionparams(rn)

    species_ids = Catalyst.speciesmap(rn)
    parameter_ids = Catalyst.paramsmap(rn)
    stoich_matrix = Catalyst.netstoichmat(rn)

    jump_rate_laws = [Catalyst.jumpratelaw(eq) for eq in Catalyst.get_eqs(rn)]

    ## Convert Catalyst jump rate laws into 
    @variables t x[1:species_count] p[1:parameter_count]

    species_sub_rules = [s => x[species_ids[s]] for s in Catalyst.species(rn)]
    parameter_sub_rules = [θ => p[parameter_ids[θ]] for θ in Catalyst.parameters(rn)]

    converted_propensities = Propensity[]
    for law in jump_rate_laws
        law1 = substitute(law, [species_sub_rules; parameter_sub_rules])
        istv = Symbolics.jacobian_sparsity([law1], [t])[]
        if istv
            b = build_function(law1, t, x, p)
        else
            b = build_function(law1, x, p)
        end
        push!(converted_propensities, propensity(eval(b)))
    end

    if typeof(parameter_values) <: AbstractVector{<:Real}
        parvec = copy(parameter_values)
    else
        pmap = Catalyst.paramsmap(rn)
        parvec = zeros(parameter_count)
        for p in parameter_values 
            parvec[ pmap[p[1]] ] = p[2]
        end
    end

    return CmeModel(stoich_matrix, converted_propensities, parvec)
end