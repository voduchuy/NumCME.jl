export propensitygrad_sparsity_pattern

function propensitygrad_sparsity_pattern(speciescount, parametercount, propensities, parameters)    
    exprs = ""
    for i in 1:length(propensities)
        if istimevarying(propensities[i])
            propensity_call = "y[$i] = propensities[$i](x[1], x[2:$speciescount+1], x[$speciescount+2:end])"                        
        else 
            propensity_call = "y[$i] = propensities[$i](x[2:$speciescount+1], x[$speciescount+2:end])"                        
        end
        exprs *= ";$propensity_call"
    end
    exprs = exprs[2:end]
    exprs = Meta.parse(exprs)
    F = eval( quote
        function prop(y, x, propensities)
            $(exprs.args...)
        end
        prop
    end)
    input = [[1.0]; ones(Int32, speciescount); parameters]
    output = zeros(length(propensities))
    sparsity_pattern = jacobian_sparsity(F, output, input, propensities, verbose=false)
    return sparsity_pattern[:, speciescount+2:end]
end