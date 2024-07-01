# Vectorized root finding

using Revise
using IntervalArithmetic
using StaticArrays
using LinearAlgebra
using ForwardDiff

# setrounding(Interval, :accurate)

IntervalArithmetic.configure!(directed_rounding=:fast, powers=:fast)

include("forwarddiff_intervalbox.jl")
include("krawczyk.jl")

# include("forwarddiff_gpu.jl")

# using CUDA
# CUDA.allowscalar(false)

function bisect_array(X)
    if length(X) == 0
        return X
    end
    
    bisected = bisect.(X)

    return vcat(first.(bisected), last.(bisected))
end



function vecroots_bisection(f, Xs::AbstractVector, numsteps=10)

    intervals = Xs

    for i in 1:numsteps

        possible_roots = contains_zero.(f.(intervals))
        to_bisect = intervals[possible_roots]

        intervals = bisect_array(to_bisect)

        @show length(intervals)
    end

    return intervals

end



function zero_test(f, Xs)
    if isempty(Xs)
        return Xs
    end

    possible_roots = contains_zero.(f.(Xs))

    if isempty(possible_roots)
        return similar(Xs, 0)
    end

    return Xs[possible_roots]
end

function krawczyk(f, Xs)
    KXs = K.(f, Xs)

    where_roots = KXs .⊂ Xs  # indices where Krawczyk proves there are roots
    no_roots = isempty.(KXs .∩ Xs)   # where Krawczyk proves there are no roots


    found_roots = Xs[where_roots]

    Xs = Xs[.!(where_roots .| no_roots)]

    return found_roots, Xs

end

"Calculate roots with given number of steps"
function simple_vecroots(f, Xs::AbstractVector, numsteps=10)

    # Xs are the list of boxes to be checked, where they may still be un-found roots

    found_roots = similar(Xs, 0)

    for i in 1:numsteps

        Xs = zero_test(f, Xs)
        isempty(Xs) && break

        new_found_roots, Xs = krawczyk(f, Xs)

        if !isempty(new_found_roots)
            if !isempty(found_roots)
                found_roots = vcat(found_roots, new_found_roots)
            else
                found_roots = new_found_roots
            end
        end

        Xs = bisect_array(Xs)

        @show length(Xs)
    end

    return found_roots, Xs

end

# Standard queue 
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






"`completed` says whether the algorithm completed while remaining under the maximum length"
function vecroots(f, Xs::AbstractVector, ϵ=1e-3, max_length=10^6; debug=true)

    # Xs are the list of boxes to be checked, where they may still be un-found roots

    found_roots = similar(Xs, 0)
    completed = true

    while !isempty(Xs) && maximum(diam.(Xs)) > ϵ

        Xs = zero_test(f, Xs)

        
        if isempty(Xs)
            completed = true
            break
        end

        # @show Xs

        new_found_roots, Xs = krawczyk(f, Xs)

        if !isempty(new_found_roots)
            if !isempty(found_roots)
                found_roots = vcat(found_roots, new_found_roots)
            else
                found_roots = new_found_roots
            end
        end


        if length(Xs) > max_length
            completed = false
            break
        end

        if length(Xs) == 0
            completed = true
            break
        end

        Xs = bisect_array(Xs)

        if debug
            @show length(Xs)
        end
    end


    if !isempty(found_roots)
        # refine
        for i in 1:10
            found_roots .= K.(f, found_roots)
            #found_roots = new_found_roots
        end
    end

    return completed, found_roots, Xs

end

# simple_vecroots(f, X, numsteps=10) = simple_vecroots(f, [X], numsteps)


# simple_vecroots(x->SVector(x[1]^2 - 2), IntervalBox(-10..10))
#
# simple_vecroots(x->SVector(x[1]^2 - 2), IntervalBox(-10..10), 2.0)

#= We need a "meta branch and bound" where we adaptively realise if a given region of space
is too complicated to solve to a given tolerance.

Branch and bound on top of branch and bound

=#

function meta_branch_and_bound(f, Xs::Vector, ϵ=0.01, max_length=10^6; cuda=false, debug=true)

    found_roots = eltype(Xs)[]
    unknown = eltype(Xs)[]

    working = copy(Xs)

    while !isempty(working)
        X = popfirst!(working)

        if diam(X) < ϵ
            push!(unknown, X)
        end

        if cuda

            completed, new_found_roots, new_unknown = vecroots(f, cu([X]), ϵ, max_length, debug=debug)

        else
            completed, new_found_roots, new_unknown = vecroots(f, [X], ϵ, max_length, debug=debug)
        end

        if completed
            if !isempty(new_found_roots)
                append!(found_roots, new_found_roots)
            end

            if !isempty(new_unknown)
                append!(unknown, new_unknown)
            end

        else
            Ys = [X]
            for i in 1:5
                Ys = bisect_array(Ys)
            end

            append!(working, Ys)
        end

    end

    return found_roots, unknown
end
