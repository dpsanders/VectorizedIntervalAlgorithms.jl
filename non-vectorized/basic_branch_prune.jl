function basic_branch_prune(f, X, ϵ=1e-3)   # no types!
    working = [X]
    results = typeof(X)[ ]

    while !isempty(working)
        X = popfirst!(working)  # first element: queue

        if diam(X) < ϵ 
            push!(results, X)
            continue
        end

        if all(0 .∈ f(X))  # vectorized
            push!(working, bisect(X)...)
        end
    end

    return results
end


