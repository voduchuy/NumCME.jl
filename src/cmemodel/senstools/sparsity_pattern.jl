export propensitygrad_sparsity_pattern

function propensitygrad_sparsity_pattern(speciescount, parametercount, propensities)    
    exprs = ""
    for i in 1:length(propensities)
        propensity_call = "y[$i] = propensities[$i](x[1], x[2:$speciescount+1], x[$speciescount+2:end])"                        
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
    input = [[1.0]; ones(Int32, speciescount); rand(parametercount)]
    output = zeros(length(propensities))
    sparsity_pattern = jacobian_sparsity(F, output, input, propensities)
    return sparsity_pattern[:, speciescount+2:end]
end