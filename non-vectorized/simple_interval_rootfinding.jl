

using Revise
using IntervalArithmetic
using StaticArrays
using LinearAlgebra
using ForwardDiff

# setrounding(Interval, :accurate)

IntervalArithmetic.configure!(directed_rounding=:fast, powers=:fast)

include("forwarddiff_intervalbox.jl")
include("krawczyk.jl")


function simple_roots(f, Xs::AbstractVector, tol=1e-3)

    # Xs are the list of boxes to be checked, where they may still be un-found roots

    found_roots = similar(Xs, 0)

    working = Xs

    T = typeof( (Xs[1], :x) )
    results = T[]

    while !(isempty(working))

        X = popfirst!(working)

        if !(contains_zero(f(X)))
            continue 
        end

        Y = K(f, X)

        if Y ⊂ X
            push!(results, (X, :unique))
            continue 
        
        else
            X = Y ∩ X
        end
        
        if diam(X) < tol
            push!(results, (X, :unknown))
            continue 
        end

        push!(working, bisect(X)...)

    end

    return results
end